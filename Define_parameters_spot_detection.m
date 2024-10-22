function Parameters = Define_parameters_spot_detection(Parameters)

prompt = {  'Use GPU?', ...
            'Substacks?', ...
            'Stack_min', ...
            'Stack_max', ...
            'Background removal?', ...
            'Background: sigma', ...
            'Intensity adjustement?', ...
            'Tolerance', ...
            'Spot detection method', ...
            'Sigma - small', ...
            'Sigma - max', ...
            'N scales', ...
            'Quantile_parameter', ...
            'h_meta_parameter', ...
            'S', ...
            'N_pixel_sample', ...
            'Perform spatial test?', ...
            };
dlgtitle = 'Parameters for seq-FISH spot detection';
dims = [1 35];

if isfield(Parameters,'use_GPU')
    definput = {
        return_boolean_text(Parameters.use_GPU), ...
        return_boolean_text(Parameters.Substack), ...
        num2str(Parameters.Stack_min), ...
        num2str(Parameters.Stack_max), ...
        return_boolean_text(Parameters.perform_background_removal), ...
        num2str(Parameters.background_sigma_parameter), ...
        return_boolean_text(Parameters.perform_intensity_adjustment), ...
        num2str(Parameters.tolerance), ...
        Parameters.Spot_detection_method, ...
        num2str(Parameters.sigma_small), ...
        num2str(Parameters.sigma_max), ...
        num2str(Parameters.N_scales), ...
        num2str(Parameters.Quantile_parameter), ...
        num2str(Parameters.h_meta_parameter), ...
        num2str(Parameters.S), ...
        num2str(Parameters.N_points_sample), ...
        return_boolean_text(Parameters.perform_spatial_statistic_test)
        };
else
    definput = {'true','false','1','25','true','80','true','0','Multiscale','1','2','3','0.00001','1','10','40000','true'};
end

answer = inputdlg(prompt,dlgtitle,dims,definput);

if ~ isempty(answer)
    
    % Read parameters
    Parameters.use_GPU = strcmp(answer{1},'true');
    Parameters.Substack = strcmp(answer{2},'true'); 
    Parameters.Stack_min = str2num(answer{3}); 
    Parameters.Stack_max = str2num(answer{4});
    Parameters.perform_background_removal = strcmp(answer{5},'true'); 
    Parameters.background_sigma_parameter = str2num(answer{6}); 
    Parameters.perform_intensity_adjustment = strcmp(answer{7},'true'); 
    Parameters.tolerance = str2num(answer{8});
	Parameters.Spot_detection_method = answer{9}; 
    Parameters.sigma_small = str2num(answer{10}); 
    Parameters.sigma_max = str2num(answer{11});
    Parameters.N_scales = str2num(answer{12}); 
    Parameters.Quantile_parameter = str2num(answer{13}); 
    Parameters.h_meta_parameter = str2num(answer{14});
    Parameters.S = str2num(answer{15}); 
    Parameters.N_points_sample = str2num(answer{16}); 
    Parameters.perform_spatial_statistic_test = strcmp(answer{17},'true'); 
else
    disp('Parameter definition canceled')
    Parameters = {};
end

