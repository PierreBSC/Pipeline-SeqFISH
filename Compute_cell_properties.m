function [Cell_properties_table] = Compute_cell_properties(Analysis_result,Parameters)
%Function that will produce basic QCs on the cells
%Provides information about cell size, DAPI staining and total spot


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


%Getting RNA and 
[Analysis_result, DAPI_staining] = Compute_cell_fluorescence(Analysis_result,Parameters,'DAPI');
[Analysis_result,Concatenated_RNA_data] = Compute_cell_RNA_expression(Analysis_result,Parameters);

%Total RNA count
Total_RNA_counts = sum(table2array(Concatenated_RNA_data),2);

%Size of the each cell
Cell_size = [];

for P=1:Parameters.N_position
    Cell_segmentation_temp = (Analysis_result.Spot_based_segmentation{P});
    Cell_size_temp = tabulate(Cell_segmentation_temp(:,3));
    Cell_size_temp = Cell_size_temp(:,2);
    Cell_size = [Cell_size;Cell_size_temp];
end

%Localisation of each cell
Cell_localisation = [];

for P=1:Parameters.N_position
    Cell_segmentation_temp = (Analysis_result.Spot_based_segmentation{P});
    Cell_localisation_temp = splitapply(@median,Cell_segmentation_temp(:,1:2),Cell_segmentation_temp(:,3));
    Cell_localisation = [Cell_localisation;Cell_localisation_temp];
end


%Image Position of each cell
Cell_position = [];

for P=1:Parameters.N_position
    Cell_segmentation_temp = (Analysis_result.Spot_based_segmentation{P});
    N_cell_temp = max(Cell_segmentation_temp(:,3));
    Cell_position = [Cell_position; repmat(P,N_cell_temp,1)];
end


Cell_properties_table = table(Cell_position,Cell_localisation,Cell_size,DAPI_staining.Concatenated_Fluo_data,Total_RNA_counts);

Cell_properties_table.Properties.RowNames =DAPI_staining.Properties.RowNames;

Cell_properties_table.Properties.VariableNames{1} = 'Position';
Cell_properties_table.Properties.VariableNames{4} = 'DAPI_intensity';
Cell_properties_table.Properties.VariableNames{5} = 'Total_RNA_molecules';

writetable(Cell_properties_table,strcat(Parameters.Output_directory,"/Cell_properties_table.txt"),'delimiter','\t','WriteRowNames',true);

end

