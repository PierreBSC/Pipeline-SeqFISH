function Parameters = Define_parameters_cell_segmentation(Parameters)
% Function that provides graphical interface to give the parameters for cell
% segmentation
prompt = {  'Min spots per cell', ...
            'N_neighbors',...
            'Number of cells ?',...
            'T (Diffusion)',...
            'Minimal similarity',...
            'Cell spread',...
            'Use filtered spots ?'
            };
dlgtitle = 'Parameters for seq-FISH cell segmentation';
dims = [1 35];

if isfield(Parameters,'Minimal_spot_cells')
    definput = {
         num2str(Parameters.Minimal_spot_cells), ...
         num2str(Parameters.N_neighbors), ...
         num2str(Parameters.N_comp), ...
         num2str(Parameters.T), ...
         num2str(Parameters.Graph_Threshold_overlap), ...
         num2str(Parameters.Bandwidth_parameter), ...
         return_boolean_text(Parameters.Perform_on_filtered_spots)
     };
else
    definput = {'20','10','30','30','0.3','15','true'};
end

answer = inputdlg(prompt,dlgtitle,dims,definput);

if ~ isempty(answer)
    % Read parameters
    Parameters.Minimal_spot_cells = str2num(answer{1}) ;
    Parameters.N_neighbors = str2num(answer{2}) ;
    Parameters.N_comp = str2num(answer{3}) ;
    Parameters.T = str2num(answer{4}) ;
    Parameters.Graph_Threshold_overlap = str2num(answer{5}) ;
    Parameters.Bandwidth_parameter = str2num(answer{6});
    Parameters.Perform_on_filtered_spots = strcmp(answer{7},'true');
    
else
    disp('Parameter definition canceled')
    Parameters = {};
end

