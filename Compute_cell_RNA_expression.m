function [Analysis_result,Concatenated_RNA_data] = Compute_cell_RNA_expression(Analysis_result,Parameters,Export_filtered_spots)
%%%Function to compute the mean intensity of each IF signal of all cells 
%%%Cell segmentation need to be performed before running this analysis 

if nargin < 3
    Export_filtered_spots = true;
end

if Parameters.Perform_on_filtered_spots
    if ~isfield(Analysis_result,'Spot_analysis_filtered')
        Parameters.Perform_on_filtered_spots =false;
    end
end

%Checking segmentation has already been performed 
if ~isfield(Analysis_result,'Spot_based_segmentation')
    disp("Perform cell segmentation step before running this function")
    return
end

%If no output directory has already been provided : provide it !

if ~isfield(Parameters,'Output_directory')
    disp("Please provide the path to the Output directory")
    Output_directory = uigetdir('title','Provide the path to the output directory');
    Parameters.Output_directory = Output_directory;
end



%Extracting spot array

if Export_filtered_spots
    Spot_array = Analysis_result.Spot_analysis_filtered;
end

if ~Export_filtered_spots
    Spot_array = Analysis_result.Spot_analysis_raw;
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

Analysis_result.RNA_analysis = cell(Parameters.N_position,1);
N_RNA_channel = sum(sum(ismember(Parameters.Matrix_design{:,:},["RNA"])));
Concatenated_RNA_data = zeros(0,N_RNA_channel);

for P=1:Parameters.N_position
    N_cells_temp = max(Analysis_result.Spot_based_segmentation{P});
    N_cells_temp = N_cells_temp(3);
    Analysis_result.Fluo_analysis{P} = zeros(N_cells_temp,N_RNA_channel);
    
end

%%Second step : creating the binary connected component corresponding for
%%each cell across each position

Cell_list = cell(Parameters.N_position,1);

for P=1:Parameters.N_position
    
    Cell_boundaries_temp = Analysis_result.Spot_based_segmentation{P,1};
    Cell_boundaries_temp(Cell_boundaries_temp==0)=1;
    N_cells_temp = max(Cell_boundaries_temp(:,3));
    M = zeros(X_size,Y_size);
    
    for i =1:N_cells_temp
        Pixels_of_interest = Cell_boundaries_temp(:,3)==i;
        M(sub2ind([size(M,1) size(M,2)],Cell_boundaries_temp(Pixels_of_interest,1),Cell_boundaries_temp(Pixels_of_interest,2))) = i;
        
    end
    Cell_list{P} = M;

end


%Across all position : extracting the cell segmentations and computing the
%number of each RNA molecule in it.

for P = 1:Parameters.N_position
        
    Temp_RNA_matrix = [];
    N_cells_temp = max(Analysis_result.Spot_based_segmentation{P});
    N_cells_temp = N_cells_temp(3);
    
    for R = 1:Parameters.N_round
    
    %If no RNA channel at this round : go to next round...
    if sum(sum(ismember(Parameters.Matrix_design{R,:},["RNA"])))==0
        continue
    end
    
    %Which are the RNA channels ?
    RNA_channel = find(ismember(Parameters.Matrix_design{R,:},["RNA"]));
    
    M =  Cell_list{P};
    
    for k = 1:size(RNA_channel,2)
        temp_spot = Spot_array{R,RNA_channel(k),P};
        temp_spot = round(temp_spot(:,1:2));
        Assigned_cell_temp =  M(sub2ind([size(M,1) size(M,2)],temp_spot(:,1),temp_spot(:,2)));
        Assigned_cell_temp = categorical(Assigned_cell_temp,1:N_cells_temp);
        RNA_count = (tabulate(Assigned_cell_temp));
        RNA_count = cell2mat(RNA_count(:,2));
        Temp_RNA_matrix = [Temp_RNA_matrix, RNA_count];
        
    end
    
    
   
    end
    Analysis_result.RNA_analysis{P} = Temp_RNA_matrix;
    Concatenated_RNA_data = [Concatenated_RNA_data ; Temp_RNA_matrix];
end

%Very last step : giving names to the cells 
%Will be with the structure Cell_##_P_## with P corresponding to the
%position.

Cell_name = [];

for P=1:Parameters.N_position
    
    N_cells_temp = max(Analysis_result.Spot_based_segmentation{P});
    N_cells_temp = N_cells_temp(3);
    temp_name = [];
    for i=1:N_cells_temp
        n = strcat("Cell_",num2str(i),"_P_",num2str(P));
        temp_name = [temp_name n];
    end
    Cell_name = [Cell_name  temp_name];
end

Concatenated_RNA_data = array2table(Concatenated_RNA_data);
Concatenated_RNA_data.Properties.RowNames = Cell_name;

Gene_names = [];

for k = 1:size(Concatenated_RNA_data,2)
    x = strcat("Gene_",num2str(k));
    Gene_names = [Gene_names, x];
end

Concatenated_RNA_data.Properties.VariableNames = Gene_names;


writetable(Concatenated_RNA_data,strcat(Parameters.Output_directory,"/RNA_expression_table.txt"),'delimiter','\t','WriteRowNames',true);

end


    