function [Analysis_result] = Spot_assignment(Analysis_result,Parameters,Perform_on_filtered_spots)
%Assign spots to cells 

Probability_threshold = 0.5;
Sigma = 10;

%By default : always perform the analysis on filtered spots, increase the
%robustness !

if nargin < 3
    Perform_on_filtered_spots = true;
end


%Checking that the spot detection and segmentation have already been
%performed : if not stop the function 

if ~isfield(Analysis_result,'Spot_analysis_raw')
     disp("Perform spot detection step before running this function")
    return
end

if ~isfield(Analysis_result,'Spot_based_segmentation')
    disp("Perform spot based segmentation step before running this function")
    return
end

if Perform_on_filtered_spots
    if ~isfield(Analysis_result,'Spot_analysis_filtered')
        Perform_on_filtered_spots =false;
    end
end


%Now we can perform cell segmentation by itself
%Creating a new field in the Analysis_result object


Analysis_result.Cell_expression = cell(Parameters.N_position,1);
Analysis_result.Cell_position = cell(Parameters.N_position,1);


%Extracting spot array

if Perform_on_filtered_spots
    Spot_array = Analysis_result.Spot_analysis_filtered;
end

if ~Perform_on_filtered_spots
    Spot_array = Analysis_result.Spot_analysis_raw;
end



for P = 1:Parameters.N_position
    
    Temp_meta_clusters = Analysis_result.Spot_based_segmentation{P,1}; %Extracting the clustering for this position
    Temp_spot_list = Spot_array(:,:,P); %Extracting the spots from the position of interest
    Temp_spot_array = [];
    
    %Aggregating all the spots from the same position
    for R =1:Parameters.N_round
        for k=1:Parameters.N_channel
            
            X = Temp_spot_list{R,k};
            
            if size(X,1)>0
                Spot_Round = repmat(R,size(X,1),1);
                Spot_Channel = repmat(k,size(X,1),1);
                X = [X,Spot_Round,Spot_Channel];
                Temp_spot_array = [Temp_spot_array ; X];
            end       
            
        end
    end
    
    N_spots_temp = size(Temp_spot_array,1);
    Assigned_cell = [];
    X = [];
    for N = 1:N_spots_temp
        
        temp_position = Temp_spot_array(N,1:2);
        Similarity = pdist2(Temp_meta_clusters(:,1:2),temp_position);
        Similarity = exp(-Similarity.^2/Sigma^2);
        Similarity = splitapply(@mean,Similarity,Temp_meta_clusters(:,4));
        Similarity = Similarity./sum(Similarity);
        [Max_similarity,Assigned_cell_temp] = max(Similarity);
        X= [X,Max_similarity];
        if Max_similarity<Probability_threshold
            Assigned_cell_temp = 0;
        end
        Assigned_cell = [Assigned_cell Assigned_cell_temp];
        

    end
    
end

end

