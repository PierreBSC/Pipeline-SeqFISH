function [Processed_image]=Pre_processing(Concat_Image,use_GPU,Substack,Stack_min,Stack_max,perform_background_removal,background_sigma_parameter,perform_intensity_adjustment,tolerance)
%Pre_processing : clean and pre-process the loaded images. 
%This step is extremely important for the spot detection so be extremely
%carefull with the settings used !
%The use_GPU option allows to perform the analysis on the computer GPU,
%hence accelerating the speed greatly


if nargin < 2
    use_GPU = true;
    Substack = false;
    perform_background_removal = true;
    background_sigma_parameter = 30;
    perform_intensity_adjustment = true;
    tolerance = 0;
end



%%Removing out of focus stacks
if Substack
    Concat_Image = Concat_Image(:,:,Stack_min:Stack_max,:);
end

numstack = size(Concat_Image,3);
l = size(Concat_Image,4);

%%If needed sending the Image data to the GPU memory
if use_GPU
    if gpuDeviceCount > 0
        Concat_Image = gpuArray(Concat_Image);
        disp("Using GPU computing for pre-processing...");
    end
end


%%Performing background removal using Gaussian smoothing
if perform_background_removal
    disp("Perfoming background removal");
    for k = 1:l 
        for i = 1:numstack
        bg_image = imgaussfilt(Concat_Image(:,:,i,k),background_sigma_parameter);
        Concat_Image(:,:,i,k) = Concat_Image(:,:,i,k) - bg_image;
        end
    end
end

Concat_Image(Concat_Image<0) = 0;


%%Performing intensity adjustment 

if perform_intensity_adjustment
    disp("Perfoming contrast/intensity adjustment.");
    for k = 1:l 
        %Geting upper and lower bound for the new intensity scale
        temp_limit = Concat_Image(:,:,:,k);
        temp_limit = stretchlim(temp_limit(:),tolerance);
        %Applying the new scale to the picture, stack by stack
        for i = 1:numstack
                    Concat_Image(:,:,i,k)=imadjust(Concat_Image(:,:,i,k),temp_limit);
        end
    end
end




%%Geting the image data back to the main memory
if use_GPU
    if gpuDeviceCount > 0
        Concat_Image = gather(Concat_Image);
        disp("Collecting data from GPU memory.")
    end
end

Processed_image=single(Concat_Image);
disp("Pre-processing finished !");



end

