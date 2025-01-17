function [] = Segmentation_visualisation(Analysis_result,Parameters,P,R,Channel)
%Visualisation function to view resulst of analysis segmentation


Round_directory = strcat(Parameters.Image_directory,"/Round_",string(R),"/");
Round_directory = char(Round_directory);
Position_directory = strcat(Round_directory,"/Position_",string(P),"/");
Position_directory = char(Position_directory);

RNA_data=LoadImage(Position_directory,true,Channel); 
RNA_data = Pre_processing(RNA_data,Parameters.use_GPU,Parameters.Substack,Parameters.Stack_min,Parameters.Stack_max,Parameters.perform_background_removal,Parameters.background_sigma_parameter,Parameters.perform_intensity_adjustment,Parameters.tolerance); 

focus_score = [];
for k = 1:size(RNA_data,3)
    focus_score = [focus_score fmeasure(RNA_data(:,:,k),'HELM')];
end
[~ , best_stack] = max(focus_score);
RNA_data = RNA_data(:,:,best_stack);


if Parameters.Perform_on_filtered_spots
    Temp_spots_list = Analysis_result.Spot_analysis_filtered(:,:,P);
end

if ~Parameters.Perform_on_filtered_spots
    Temp_spots_list = Analysis_result.Spot_analysis_raw(:,:,P);
end

Temp_spots_list = Temp_spots_list(:);

%Let's plot the spots from the different channels 

Aggregated_points = [];

for i=1:(size(Temp_spots_list,1))
    x = Temp_spots_list{i};
    if size(x,1)>0
        x = x(:,1:2);
        x = [x,repmat(i,size(x,1),1)];
        Aggregated_points = [Aggregated_points;x];
    end
end


figure, imshow(imadjust(RNA_data,stretchlim(RNA_data,0.001)))
hold on
gscatter(Aggregated_points(:,1),Aggregated_points(:,2),Aggregated_points(:,3))

%Now the enveloppes
Selected_segmentation = Analysis_result.Spot_based_segmentation{P};
N_cells = size(unique(Selected_segmentation(:,3)),1);

for k=1:N_cells
    X = Selected_segmentation(Selected_segmentation(:,3)==k,[1 2]);
    Final_boundary = X(boundary(X(:,1:2)),:);
    plot(Final_boundary(:,2),Final_boundary(:,1), 'LineWidth',3)
end

end
