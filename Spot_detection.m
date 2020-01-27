function Analysis_result=Spot_detection(Analysis_result,Parameters)
%Generic function that performs spot detection
%Requires an Analysis_result and Parameters objects to already exist

disp(' >>> Start spot detection')

if Parameters.use_GPU
    if gpuDeviceCount > 0
        d = gpuDevice;
        fprintf('\nWill use GPU computing for some analysis steps: %s.\n',d.Name)
    end
end


%First step : create a new element inside the Analysis_result to store the
%results of spot detection analysis

Analysis_result.Spot_analysis_raw = cell(Parameters.N_round,Parameters.N_channel,Parameters.N_position);

if Parameters.perform_spatial_statistic_test
    Analysis_result.Spot_analysis_filtered = cell(Parameters.N_round,Parameters.N_channel,Parameters.N_position);
end

for R = 1:Parameters.N_round
    
    %If no RNA channel at this round : go to next round...
    if sum(sum(Parameters.Matrix_design{R,:}=="RNA"))==0
        continue
    end
    
    %Defining the Round directory
    Round_directory = strcat(Parameters.Image_directory,'/Round_',string(R),'/');
    Round_directory = char(Round_directory);
    
    for P = 1:Parameters.N_position
        
        fprintf(strcat("\n >>> Processing: Round ",num2str(R)," Position ",num2str(P),"\n"))
        
        % Defining the Position subdirectory
        Position_directory = strcat(Round_directory,'/Position_',string(P),'/');
        Position_directory = char(Position_directory);
        
        % Which are the RNA channels ?
        RNA_channel = find(Parameters.Matrix_design{R,:}=="RNA");
        
        % Loading only RNA data
        RNA_data=LoadImage(Position_directory,true,RNA_channel);
        
        %Cleaning of the data
        RNA_data = Pre_processing(RNA_data,Parameters.use_GPU,Parameters.Substack,Parameters.Stack_min,Parameters.Stack_max,Parameters.perform_background_removal,Parameters.background_sigma_parameter,Parameters.perform_intensity_adjustment,Parameters.tolerance);     %Detection through the Multiscale or H-dome method
        
        % Perform spot detection with selected method
        switch Parameters.Spot_detection_method
            
            case "Multiscale"
                disp('Performing spot detection using Multiscale method')
                Unfiltered_spots = Multiscale_filter(RNA_data,Parameters.use_GPU,Parameters.sigma_small,Parameters.sigma_max,Parameters.N_scales,Parameters.T_offset);
                
            case "H-dome"
                disp('Performing spot detection using HD method')
                Unfiltered_spots = H_dome_filter(RNA_data,Parameters.use_GPU,Parameters.sigma_small,Parameters.sigma_max,Parameters.h_meta_parameter,Parameters.S,Parameters.N_points_sample);
                
            otherwise
                disp('No known spot detection method provided. Using Multiscale method  by default.')
                Unfiltered_spots = Multiscale_filter(RNA_data,Parameters.use_GPU,Parameters.sigma_small,Parameters.sigma_max,Parameters.N_scales,Parameters.T_offset);
        end
       
        Analysis_result.Spot_analysis_raw(R,RNA_channel,P) = Unfiltered_spots; %Exporting the raw spots
        
        % If specified, a spatial statistical test is performed
        % to filter out signal of non-specificly  bound probes
        if Parameters.perform_spatial_statistic_test
            temporary_directory = '~/seqFISH_Temporary_file/';
            [Intensity_threshold, Spot_intensity] = Spot_parameter_optimisation(RNA_data,Unfiltered_spots,Parameters.use_GPU,temporary_directory,strcat(pwd,'/Scan_test.R'),Parameters.sigma_small); %Computing the value intensity threshold
            Filtered_spots = Spot_filtering(Unfiltered_spots,Spot_intensity,Intensity_threshold); %%Filtering of the spots
            Analysis_result.Spot_analysis_filtered(R,RNA_channel,P) = Filtered_spots; %%
        end
        
        %clear RNA_data
        
    end
end

disp('End of the spot detection analysis')



