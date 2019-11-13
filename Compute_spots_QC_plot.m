function [] = Compute_spots_QC_plot(Analysis_result,Parameters,temporary_directory)
%Function creating the QC plots corresponding to the spot analysis 
%   Compute multiple statistics of the spots using R package spatstat
%   functions and graphical abilities...

if nargin < 3
    temporary_directory = '~/seqFISH_Temporary_file/';
end

%If the temporary file does not exist : create it !
if ~exist(temporary_directory,'dir')
    mkdir(temporary_directory)
end

%Remove all files from the directory
delete(strcat(temporary_directory,'/*'))

%Concatenating the spot information across all Rounds, Position and
%Channels
%Information collected : spot localisation, spot intensity, filtered or not

Concatenated_spots = [];
Estimated_noise_matrix = zeros(size(Analysis_result.Spot_analysis_raw));


for R = 1:Parameters.N_round
    
    %If no RNA channel at this round : go to next round...
    if sum(sum(Parameters.Matrix_design{R,:}=="RNA"))==0
        continue
    end
    
    %Defining the Round directory
    Round_directory = strcat(Parameters.Image_directory,'/Round_',string(R),'/');
    Round_directory = char(Round_directory);
    List_spots = [];
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

                %Computing mean intensity of the spots as well as
                %estimated background noise 
                Channel_data=LoadImage(Position_directory,true,k); 
                Channel_data = Pre_processing(Channel_data,Parameters.use_GPU,Parameters.Substack,Parameters.Stack_min,Parameters.Stack_max,Parameters.perform_background_removal,Parameters.background_sigma_parameter,Parameters.perform_intensity_adjustment,Parameters.tolerance); 
                
                Estimated_noise = std(Channel_data(:));
                Estimated_noise_matrix(R,k,P) = Estimated_noise;
                Smoothed_image = imgaussfilt3(Channel_data,Parameters.sigma_small);
                Estimated_intensity = Smoothed_image(sub2ind(size(Smoothed_image),round(X_raw(:,1)),round(X_raw(:,2)),round(X_raw(:,3))));
                 
                %Final table : position (X,Y,Z),
                X = [X_raw Estimated_intensity Is_kept repmat(R,size(X_raw,1),1) repmat(P,size(X_raw,1),1) repmat(k,size(X_raw,1),1)];
                List_spots = [List_spots; X];

            end
            
        end
    end
end

%%Changing the format of the Noise matrix so it can be loaded by R

if size(Estimated_noise_matrix,3)==1
[x_index y_index z_index] = ind2sub([size(Estimated_noise_matrix) 1],1:size(Estimated_noise_matrix(:),1));
end

if size(Estimated_noise_matrix,3)>1
[x_index y_index z_index] = ind2sub(size(Estimated_noise_matrix),1:size(Estimated_noise_matrix(:),1));
end

Estimated_noise_table = [Estimated_noise_matrix(:) , x_index.' , y_index.', z_index.'];
%Exporting the final table in the temporary directory

writetable(array2table(List_spots),strcat(temporary_directory,'/Concatenated_spot_information.txt'),'Delimiter','\t');
writetable(array2table(Estimated_noise_table),strcat(temporary_directory,'/Noise_value.txt'),'Delimiter','\t');

%Exporting basic information about the experiment settings (number of
%rounds, positions and channels)

writetable(array2table([Parameters.N_round Parameters.N_channel Parameters.N_position]),strcat(temporary_directory,'/General_experiment_information.txt'),'Delimiter','\t');





end

