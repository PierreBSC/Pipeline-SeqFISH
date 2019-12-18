function [Parameters] = Compute_spots_QC_plot(Analysis_result,Parameters,temporary_directory,noise_box)
%Function creating the QC plots corresponding to the spot analysis 
%   Compute multiple statistics of the spots using R package spatstat
%   functions and graphical abilities...

if nargin < 3
    temporary_directory = '~/seqFISH_Temporary_file/';
    noise_box = [3 3]; %What is the size of the patch used to compute local std ? Here corresponds to a 3+1+3 * 3+1+3 box...
end

QC_spots_script = strcat(pwd,'/QC_spots_script.R');

%If the temporary file does not exist : create it !
if ~exist(temporary_directory,'dir')
    mkdir(temporary_directory)
end

%Remove all files from the directory
delete(strcat(temporary_directory,'/*'))

%If no output directory has already been provided : provide it !

if ~isfield(Parameters,'Output_directory')
    disp("Please provide the path to the Output directory")
    Output_directory = uigetdir('title','Provide the path to the output directory');
    Parameters.Output_directory = Output_directory;
end


%%We need to know what is the size of the pictures : just load one
%%picture...
Round_directory = strcat(Parameters.Image_directory,"/Round_1/");
Round_directory = char(Round_directory);
Position_directory = strcat(Round_directory,"/Position_1/");
Position_directory = char(Position_directory);
Example_data=LoadImage(Position_directory,true,1); 
X_size = size(Example_data,1);
Y_size = size(Example_data,2);


%Concatenating the spot information across all Rounds, Position and
%Channels
%Information collected : spot localisation, spot intensity, filtered or not
List_spots = [];



for R = 1:Parameters.N_round
    
    %If no RNA channel at this round : go to next round...
    if sum(sum(Parameters.Matrix_design{R,:}=="RNA"))==0
        continue
    end
    
    %Defining the Round directory
    Round_directory = strcat(Parameters.Image_directory,'/Round_',string(R),'/');
    Round_directory = char(Round_directory);
    
    for P = 1:Parameters.N_position
    
        List_spots_temp = [];
        Position_directory = strcat(Round_directory,'/Position_',string(P),'/');
        Position_directory = char(Position_directory);

        for k = 1:size(Parameters.Matrix_design,2)
            
            if size(Analysis_result.Spot_analysis_raw{R,k,P},1) > 0
                
                X_raw = Analysis_result.Spot_analysis_raw{R,k,P};
                
                if Parameters.perform_spatial_statistic_test
                    X_filtered = Analysis_result.Spot_analysis_filtered{R,k,P};
                    Is_kept = ismember(X_raw,X_filtered,'row');
                end
                
                 if ~Parameters.perform_spatial_statistic_test
                    Is_kept = repmat(0,size(X_raw,1),1);
                end

                %Computing intensity of the spots as well as
                %estimated local background noise f
                Channel_data=LoadImage(Position_directory,true,k); 
                Channel_data = Pre_processing(Channel_data,Parameters.use_GPU,Parameters.Substack,Parameters.Stack_min,Parameters.Stack_max,Parameters.perform_background_removal,Parameters.background_sigma_parameter,Parameters.perform_intensity_adjustment,Parameters.tolerance); 

                Estimated_intensity = Channel_data(sub2ind(size(Channel_data),round(X_raw(:,1)),round(X_raw(:,2)),round(X_raw(:,3))));
                Estimated_noise = [];
                
                for i=1:size(X_raw,1)
                    Point_center = [round(X_raw(i,1)) round(X_raw(i,2)) round(X_raw(i,3))];
                    X_min = max(1,Point_center(1)-noise_box(1));
                    Y_min = max(1,Point_center(2)-noise_box(2));
                    X_max = min(X_size,Point_center(1)+noise_box(1));
                    Y_max = min(Y_size,Point_center(2)+noise_box(2));
                    [List_sub_x,List_sub_y] = meshgrid(X_min:X_max,Y_min:Y_max);
                    Channel_data_temp = Channel_data(:,:,round(X_raw(i,3)));
                    Local_values = Channel_data_temp(sub2ind(size(Channel_data_temp),List_sub_x(:),List_sub_y(:)));
                    Estimated_noise =[Estimated_noise std(Local_values)];
                end
                 
                %Final table : position (X,Y,Z),
                X = [X_raw Estimated_intensity Is_kept Estimated_noise.' repmat(R,size(X_raw,1),1) repmat(P,size(X_raw,1),1) repmat(k,size(X_raw,1),1)];
                List_spots = [List_spots; X];

            end
            
        end
    end
end
writetable(array2table(List_spots),strcat(temporary_directory,'/Concatenated_spot_information.txt'),'Delimiter','\t');

%Exporting basic information about the experiment settings (number of
%rounds, positions and channels)

writetable(array2table([Parameters.N_round Parameters.N_channel Parameters.N_position]),strcat(temporary_directory,'/General_experiment_information.txt'),'Delimiter','\t');


system(strcat("Rscript"," ",QC_spots_script," ",temporary_directory, " ", Parameters.Output_directory ));


end

