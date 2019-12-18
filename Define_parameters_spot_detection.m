function Parameters = Define_parameters_spot_detection()

prompt = {  'Name of experiment', ...
            'Number of positions', ...
            'Use GPU?', ...
            'Substacks?', ...
            'Stack_min', ...
            'Stack_max', ...
            'Background removal?', ...
            'Background: sigma', ...
            'Intensity adjustement?', ...
            'Tolerance', ...
            'Spot detection method', ...
            'Sigma - small', ...
            'N scales', ...
            'T offset', ...
            'Sigma - max', ...
            'h_meta_parameter', ...
            'S', ...
            'N_pixel_sample', ...
            'Perform spatial test?', ...
            };
dlgtitle = 'Parameters for seq-FISH spot detection';
dims = [1 35];
definput = {'seq-FISH analysis','1','true','false','1','25','true','20','true','0','Multiscale','1','10','0','8','1.5','10','40000','true'};

answer = inputdlg(prompt,dlgtitle,dims,definput);

if ~ isempty(answer)
    
    % Read parameters
    Parameters.Experiment_name = answer{1}; 
    Parameters.N_position = str2num(answer{2});
    Parameters.use_GPU = strcmp(answer{3},'true');
    Parameters.Substack = strcmp(answer{4},'true'); 
    Parameters.Stack_min = str2num(answer{5}); 
    Parameters.Stack_max = str2num(answer{6});
    Parameters.perform_background_removal = strcmp(answer{7},'true'); 
    Parameters.background_sigma_parameter = str2num(answer{8}); 
    Parameters.perform_intensity_adjustment = strcmp(answer{9},'true'); 
    Parameters.tolerance = str2num(answer{10});
	Parameters.Spot_detection_method = answer{11}; 
    Parameters.sigma_small = str2num(answer{12}); 
    Parameters.N_scales = str2num(answer{13}); 
    Parameters.T_offset = str2num(answer{14}); 
    Parameters.sigma_max = str2num(answer{15});
    Parameters.h_meta_parameter = str2num(answer{16});
    Parameters.S = str2num(answer{17}); 
    Parameters.N_points_sample = str2num(answer{18}); 
    Parameters.perform_spatial_statistic_test = strcmp(answer{19},'true'); 
    
else
    disp('Parameter definition canceled')
    Parameters = {};
end
