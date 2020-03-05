

function [Intensity_threshold Spot_intensity Autolag_score] = Spot_parameter_optimisation(Processed_Image,Unfiltered_spots,use_GPU,temporary_directory,Spatial_statistic_script,Estimated_size_cell)


%Estimation of a different intensity thresold parameter for each channel using the Kulldorf scan test 

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

for k=1:l

    u = Processed_Image(:,:,:,k);

    channel_spot = round(Unfiltered_spots{k});

    Spot_intensity{k}=table(gather(u(sub2ind(size(u),channel_spot(:,2),channel_spot(:,1),channel_spot(:,3)))));

end




for k=1:l 

    writetable((Spot_intensity{k}),strcat(temporary_directory,'Spot_intensity_channel_',num2str(k),'.txt'),'Delimiter','\t');
    writetable(table(Unfiltered_spots{k}),strcat(temporary_directory,'Spot_position_channel_',num2str(k),'.txt'),'Delimiter','\t');

end

system(strcat("Rscript"," ",Spatial_statistic_script," ",temporary_directory," ",num2str(Estimated_size_cell)));


Intensity_threshold = zeros(1,l);


for k=1:l

    table_threshold = readtable(strcat(temporary_directory,'Table_channel',num2str(k),'.txt'));

    [x ,t] = min(table_threshold.Score);

    Intensity_threshold(k) = table_threshold.Intensity(max(t,1));

end

Autolag_score = readmatrix(strcat(temporary_directory,'Table_autolag_score.txt'));


disp('Automated parameter estimation done')


end


