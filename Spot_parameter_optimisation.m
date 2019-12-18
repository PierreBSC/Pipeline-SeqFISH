

function [Intensity_threshold Spot_intensity] = Spot_parameter_optimisation(Processed_Image,Unfiltered_spots,use_GPU,temporary_directory,Spatial_statistic_script,Estimated_size_spot)


%Estimation of a different intensity thresold parameter for each channel using the Kulldorf scan test 

if nargin < 5Round_directory = strcat(Parameters.Image_directory,"/Round_",string(R),"/");
Round_directory = char(Round_directory);
Position_directory = strcat(Round_directory,"/Position_",string(P),"/");
Position_directory = char(Position_directory);

RNA_data=LoadImage(Position_directory,true,Channel); 
RNA_data = Pre_processing(RNA_data,Parameters.use_GPU,Parameters.Substack,Parameters.Stack_min,Parameters.Stack_max,Parameters.perform_background_removal,Parameters.background_sigma_parameter,Parameters.perform_intensity_adjustment,Parameters.tolerance); 
%RNA_data = mean(RNA_data,3);
focus_score = [];


    temporary_directory = '~/seqFISH_Temporary_file/';

    use_GPU = true;

    Estimated_size_spot = 1;

    Spatial_statistic_script = strcat(pwd,'/Scan_test.R');

end


%If the temporary file does not exist : create it !

if ~exist(temporary_directory,'dir')

    mkdir(temporary_directory)

end


%Remove all files from the directory

delete(strcat(temporary_directory,'/*'))


%%If needed sending the Image data to the GPU memory

if use_GPU

    if gpuDeviceCount > 0

        Processed_Image = gpuArray(Processed_Image);

    end

end




%%First step : getting spot intensity 


l = size(Processed_Image,4);

Intensity_threshold = zeros(l,1);


Spot_intensity=cell(l,1); %%One list containing the intensity of the chanel in which the spot was detected 


%%%We start by smoothing the images using a 3D Gaussian Filter 

%%%We will use the smooth intensity as the corrected mean intensity of the

%%%probe


Smoothed_image = single(zeros(size(Processed_Image)));



if use_GPU

    if gpuDeviceCount > 0

        Smoothed_image = gpuArray(Smoothed_image);

    end

end



for k=1:l

    Smoothed_image(:,:,:,k) = imgaussfilt3(Processed_Image(:,:,:,k),Estimated_size_spot);

end


%%%We then copy the values in the different variables

for k=1:l

    u = Smoothed_image(:,:,:,k);

    channel_spot = round(Unfiltered_spots{k});

    Spot_intensity{k}=table(gather(u(sub2ind(size(u),channel_spot(:,2),channel_spot(:,1),channel_spot(:,3)))));

end




for k=1:l 

    writetable((Spot_intensity{k}),strcat(temporary_directory,'Spot_intensity_channel_',num2str(k),'.txt'),'Delimiter','\t');

    writetable(table(Unfiltered_spots{k}),strcat(temporary_directory,'Spot_position_channel_',num2str(k),'.txt'),'Delimiter','\t');

end

system(strcat("Rscript"," ",Spatial_statistic_script," ",temporary_directory));


Intensity_threshold = zeros(1,l);


for k=1:l

    table_threshold = readtable(strcat(temporary_directory,'Table_channel',num2str(k),'.txt'));

    [~ ,t] = max(table_threshold.Score);

    Intensity_threshold(k) = table_threshold.Intensity(max(t,1));

end


disp('Automated parameter estimation done')


end


