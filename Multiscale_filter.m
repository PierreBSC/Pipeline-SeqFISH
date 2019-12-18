function [List_spots] = Multiscale_filter(Processed_Image,use_GPU,sigma_small,N,Threshold_adjustement)
%Spot detection performed using multiscale filtering
%The 



if nargin < 5
    use_GPU = true;
    sigma_small = 1; %% Size of the smallest gaussian filter, in pixel 
    N = 3 ; %%Number of scales used 
    Threshold_adjustement = 0.05;
    disp('All needed parameters have not been specified. Running using default parameters')
end



%%If needed sending the Image data to the GPU memory
if use_GPU
    if gpuDeviceCount > 0
        Processed_Image = gpuArray(Processed_Image);
        disp('Using GPU computing for pre-processing...');
    end
end

%%Extracting image information
nImage = size(Processed_Image,1);
mImage = size(Processed_Image,2);
numstack = size(Processed_Image,3);
l = size(Processed_Image,4);


%Creating objects
sigma_list =  sigma_small * 2.^(0:(N-1));
List_matrix = zeros([size(Processed_Image),N]);


%Computing the determinant of the Hessian matrix at different scales on
%each stack for each channel

disp("Computing Hessian Matrix...")
for k=1:l
    for i=1:numstack
        for j = 1:N
            Smoothed_data = imgaussfilt(Processed_Image(:,:,i,k),sigma_list(j));
            [gx, gy] = imgradientxy(Smoothed_data);
            [gxx, gxy] = imgradientxy(gx);
            [gyx, gyy] = imgradientxy(gy);
            Determinant = gxx.*gyy - 2*gxy;
            List_matrix(:,:,i,k,j) = gather(Determinant);
        end
    end
end

%Extracting the highest determinant obtained 


Cleaned_signal = max(List_matrix,[],5);

%Performing discretization using Triangle method

disp("Performing automated thresholding...")

Discretized_signal = zeros(size(Processed_Image));

parfor k = 1:l
    
    Cleaned_signal_temp = (Cleaned_signal(:,:,:,k));
    Cleaned_signal_temp = Cleaned_signal_temp./max(Cleaned_signal_temp(:));
    Cleaned_signal_list = Cleaned_signal_temp(:);
    
    [Histogram_count Histogram_bin] = hist(Cleaned_signal_list,100);
    Peak_value = Histogram_bin(find(Histogram_count==max(Histogram_count)));
    Histogram_count = Histogram_count(Histogram_bin >= Peak_value);
    Histogram_bin = Histogram_bin(Histogram_bin >= Peak_value);
    
    Histogram_count = Histogram_count./sum(Histogram_count);
    Histogram_line = 1+Histogram_bin.*(1/(Peak_value-1));
    Line_points = [Histogram_bin.',Histogram_line.'];
    Histo_points = [Histogram_bin.',Histogram_count.'];
    
    [X, dist_closest] = knnsearch(Line_points,Histo_points,'K',1);
    Threshold = Histogram_bin(find(dist_closest==max(dist_closest)));

    Discretized_signal(:,:,:,k) = Cleaned_signal_temp > Threshold+ Threshold_adjustement;
end


%Identifying the connected component

disp("Connected components extraction and filtering")

List_spots = cell(l,1);

parfor k=1:l   
    
    List_spots_temp = bwconncomp(Discretized_signal(:,:,:,k),6);
    Size_spots = regionprops3(List_spots_temp,'Volume');
    Size_spots = Size_spots.Volume;
    X = Cleaned_signal(:,:,:,k);
    Position_spots = regionprops3(List_spots_temp,gather(X),'WeightedCentroid');
    Position_spots = Position_spots.WeightedCentroid;
    Position_spots_filtered = Position_spots(Size_spots>5,:);
    List_spots{k}=Position_spots_filtered;
        
end


end

