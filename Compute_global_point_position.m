function [Global_analysis_results] = Compute_global_point_position(Analysis_result,Parameters)
%Compute an intensity threshold using all pictures stitched together

if ~isfield(Analysis_result,'Global_stitching')
    disp("Perform image stitching step before running this function")
    return
end

%Creating the different fields
Global_analysis_results = struct;
Global_analysis_results.Spot_analysis_raw = cell(Parameters.N_round,Parameters.N_channel);

if Parameters.perform_spatial_statistic_test
    Global_analysis_results.Spot_analysis_filtered = cell(Parameters.N_round,Parameters.N_channel);
end

%%Aggregation of the spot positions

for R = 1:Parameters.N_round
    
    for P = 1:Parameters.N_position
    disp(num2str(P))

       for k=1:Parameters.N_channel
           
        %Let's start with the unfiltered spots  
        Spot_position_raw_temp = Analysis_result.Spot_analysis_raw{R,k,P};
        
        if Parameters.perform_spatial_statistic_test
            Spot_position_filtered_temp = Analysis_result.Spot_analysis_filtered{R,k,P};
        end

        T = affine2d([1 0 0 ; 0 1 0 ; -Analysis_result.Global_stitching(P,1) ,-Analysis_result.Global_stitching(P,2) 1]);

        if size(Spot_position_raw_temp,1) > 0
            
            Spot_position_raw_temp = [T.transformPointsForward(Spot_position_raw_temp(:,1:2)),Spot_position_raw_temp(:,3)];
            
            %We need to remove duplicate points
            
            if P > 1
            X = Global_analysis_results.Spot_analysis_raw{R,k};
            [closest_point,closest_dist]= knnsearch(X(:,1:2),Spot_position_raw_temp(:,1:2),'K',1);
            N_shared_points = sum(closest_dist < Parameters.sigma_small*2);
            List_new_point_ID = find(closest_dist < Parameters.sigma_small*2);
            Spot_position_raw_temp(List_new_point_ID,:) = [];
                            
            end            

        end
        
        Spot_position_filtered_temp = Analysis_result.Spot_analysis_filtered{R,k,P};

        if size(Spot_position_filtered_temp,1) > 0
            
            Spot_position_filtered_temp = [T.transformPointsForward(Spot_position_filtered_temp(:,1:2)),Spot_position_filtered_temp(:,3)];
            %We need to remove duplicate points
            
            if P>1
            %We need to remove duplicate points
            
            X = Global_analysis_results.Spot_analysis_filtered{R,k};
            [closest_point,closest_dist]= knnsearch(X(:,1:2),Spot_position_filtered_temp(:,1:2),'K',1);
            N_shared_points = sum(closest_dist < Parameters.sigma_small*2);
            List_new_point_ID = find(closest_dist < Parameters.sigma_small*2);
            Spot_position_filtered_temp(List_new_point_ID,:) = [];

            end
        
        
        end
        
        %Saving the files
        Global_analysis_results.Spot_analysis_raw{R,k} = [Global_analysis_results.Spot_analysis_raw{R,k} ; Spot_position_raw_temp];
        Global_analysis_results.Spot_analysis_filtered{R,k} = [Global_analysis_results.Spot_analysis_filtered{R,k} ; Spot_position_filtered_temp];

       end
    end
end

%%Filtering the spots that are too close






end

