function [] = Spot_visualisation(Analysis_result,Parameters,P,R,Channel)
%Visualisation function to view resulst of spot analysis

if nargin < 3
    R = 1;
    P = 1;
    Channel = 1;
end

Round_directory = strcat(Parameters.Image_directory,"/Round_",string(R),"/");
Round_directory = char(Round_directory);
Position_directory = strcat(Round_directory,"/Position_",string(P),"/");
Position_directory = char(Position_directory);

RNA_data=LoadImage(Position_directory,true,Channel); 
RNA_data = Pre_processing(RNA_data,Parameters.use_GPU,Parameters.Substack,Parameters.Stack_min,Parameters.Stack_max,Parameters.perform_background_removal,Parameters.background_sigma_parameter,Parameters.perform_intensity_adjustment,Parameters.tolerance); 

%RNA_data = mean(RNA_data,3);
focus_score = [];
for k = 1:size(RNA_data,3)
    focus_score = [focus_score fmeasure(RNA_data(:,:,k),'ACMO')];
end
[~ , best_stack] = max(focus_score);
RNA_data = RNA_data(:,:,best_stack);

figure, imshow(imadjust(RNA_data,stretchlim(RNA_data,0.0005)))
hold on
Raw_spots = Analysis_result.Spot_analysis_raw{R,Channel,P};
scatter(Raw_spots(:,1),Raw_spots(:,2),'LineWidth',4)

if Parameters.perform_spatial_statistic_test
    Filtered_spots = Analysis_result.Spot_analysis_filtered{R,Channel,P};
    scatter(Filtered_spots(:,1),Filtered_spots(:,2),repmat(60,size(Filtered_spots,1),1),'filled')
end

end
