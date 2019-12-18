function Parameters = Define_parameters_cell_segmentation(Parameters)
%FUnction that provides graphical interface to give the parameters for cell
%segementation
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
definput = {'20','10','30','30','0.3','15','true'};

answer = inputdlg(prompt,dlgtitle,dims,definput);

if ~ isempty(answer)
    % Read parameters
    Parameters.Minimal_spot_cells = str2num(answer{1}) ;
    Parameters.N_neighbors = str2num(answer{2}) ;
    Parameters.N_comp = str2num(answer{3}) ;
    Parameters.T = str2num(answer{4}) ;
    Parameters.Graph_Threshold_overlap = str2num(answer{5}) ;
    Parameters.Bandwidth_parameter = str2num(answer{6})
    Parameters.Perform_on_filtered_spots = strcmp(answer{7},'true')
    
else
    disp('Parameter definition canceled')
    Parameters = {};
end

