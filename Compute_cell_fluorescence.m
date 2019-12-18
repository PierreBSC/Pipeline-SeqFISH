function [Analysis_result,Concatenated_Fluo_data] = Compute_cell_fluorescence(Analysis_result,Parameters,IF_or_DAPI)
%%%Function to compute the mean intensity of each IF signal of all cells 
%%%Cell segmentation needs to be performed before running this analysis 
%%% Can be used to look to either IF or DAPI stainings

if nargin < 3
    IF_or_DAPI = 'IF';
end


%%%Checking that the spots detection stage has already been performed

if ~isfield(Analysis_result,'Spot_based_segmentation')
    disp("Perform cell segmentation step before running this function")
    return
end


%%Getting size of images
%%We need to know what is the size of the pictures : just load one
%%picture...
Round_directory = strcat(Parameters.Image_directory,"/Round_1/");
Round_directory = char(Round_directory);
Position_directory = strcat(Round_directory,"/Position_1/");
Position_directory = char(Position_directory);
Example_data=LoadImage(Position_directory,true,1); 
X_size = size(Example_data,1);
Y_size = size(Example_data,2);


%%%First step : let's create the field in the Analysis_result object as well
%%%as the corresponding array 

Analysis_result.Fluo_analysis = cell(Parameters.N_position,1);
N_Fluo_channel = sum(sum(ismember(Parameters.Matrix_design{:,:},IF_or_DAPI)));
Concatenated_Fluo_data = zeros(0,N_Fluo_channel);


for P=1:Parameters.N_position
    N_cells_temp = max(Analysis_result.Spot_based_segmentation{P});
    N_cells_temp = N_cells_temp(3);
    Analysis_result.Fluo_analysis{P} = zeros(N_cells_temp,N_Fluo_channel);
    
end

%%Second step : creating the binary connected component corresponding for
%%each cell across each position

Cell_list = cell(Parameters.N_position,1);

for P=1:Parameters.N_position
    
    Cell_boundaries_temp = Analysis_result.Spot_based_segmentation{P};
    N_cells_temp = max(Cell_boundaries_temp(:,3));
    M = zeros(X_size,Y_size);
    
    for i =1:N_cells_temp
        Pixels_of_interest = Cell_boundaries_temp(:,3)==i;
        M(sub2ind([size(M,1) size(M,2)],Cell_boundaries_temp(Pixels_of_interest,1),Cell_boundaries_temp(Pixels_of_interest,2))) = i;
    end
    Cell_list{P} = M;

end


%Across all position : extracting the cell segmentations and computing the
%

for P = 1:Parameters.N_position
        
    Temp_intesity_matrix = [];
    
    for R = 1:Parameters.N_round
    
    %If no IF or DAPI channel at this round : go to next round...
    if sum(sum(ismember(Parameters.Matrix_design{R,:},IF_or_DAPI)))==0
        continue
    end

    %Defining the Round directory
    Round_directory = strcat(Parameters.Image_directory,'/Round_',string(R),'/');
    Round_directory = char(Round_directory);
    
    %Defining the Position subdirectory
    Position_directory = strcat(Round_directory,'/Position_',string(P),'/');
    Position_directory = char(Position_directory);
    
    %Which are the fluorescence channels ?
    Fluo_channel = find(ismember(Parameters.Matrix_design{R,:},IF_or_DAPI));
    
    %Loading only RNA data
    Fluo_data=LoadImage(Position_directory,true,Fluo_channel); 
    
    %No pre-processing of the IF data should be performed as it would make
    %the comparison between Positions non robust..
    % Instead : stack selection is performed for each channel individually
    
    Stacked_Fluo_data = zeros(size(Fluo_data,1),size(Fluo_data,2),size(Fluo_data,4));
    
    for k = 1:size(Fluo_data,4)
         focus_score = [];
        for i = 1:size(Fluo_data,3)
            focus_score = [focus_score fmeasure(Fluo_data(:,:,i,k),'HELM')];
        end
        [~ , best_stack] = max(focus_score);
        Stacked_Fluo_data(:,:,k) = Fluo_data(:,:,best_stack,k);
    end
    
    %%No we can extract the mean intensity of the channels
    
    for k = 1:size(Fluo_data,4)
        X = regionprops(Cell_list{P},Stacked_Fluo_data(:,:,k),'MeanIntensity');
        X = struct2array(X).';
        Temp_intesity_matrix = [Temp_intesity_matrix, X];
    end
   
    end
    Analysis_result.Fluo_analysis{P} = Temp_intesity_matrix;
    Concatenated_Fluo_data = [Concatenated_Fluo_data ; Temp_intesity_matrix];
end

%Very last step : giving names to the cells 
%Will be with the structure Cell_##_P_## with P corresponding to the
%position.

Cell_name = [];

for P=1:Parameters.N_position
    
    N_cells_temp = max(Analysis_result.Spot_based_segmentation{P,1});
    N_cells_temp = N_cells_temp(3);
    temp_name = [];
    for i=1:N_cells_temp
        n = strcat("Cell_",num2str(i),"_P_",num2str(P));
        temp_name = [temp_name n];
    end
    Cell_name = [Cell_name  temp_name];
end

Concatenated_Fluo_data = array2table(Concatenated_Fluo_data);

if size(Concatenated_Fluo_data,1) >0
    Concatenated_Fluo_data.Properties.RowNames = Cell_name;
end

if IF_or_DAPI=="IF"
    writetable(Concatenated_Fluo_data,strcat(Parameters.Output_directory,"/IF_expression_table.txt"),'delimiter','\t','WriteRowNames',true);
end


end


    