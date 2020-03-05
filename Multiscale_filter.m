function [List_spots] = Multiscale_filter(Processed_Image,use_GPU,sigma_small,sigma_max,N,Quantile_parameter)
%Spot detection performed using multiscale filtering

% Default paramter when not enough input is provided
if nargin < 5
    use_GPU = true;
    sigma_small = 1; %% Size of the smallest gaussian filter, in pixel 
    N = 4 ; %%Number of scales used 
    Quantile_parameter = 0.001;
    disp('All needed parameters have not been specified. Running using default parameters')
end

% If specified: send image data to the GPU memory
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
sigma_list =  linspace(sigma_small,sigma_max,N);
List_matrix = zeros([size(Processed_Image,1),size(Processed_Image,2),size(Processed_Image,3),size(Processed_Image,4),N]);

%Computing the determinant of the Hessian matrix at different scales on
%each stack for each channel

disp("Computing Hessian Matrix...")
for k=1:l
    for i=1:numstack
        for j = 1:N
            Smoothed_data = imgaussfilt(Processed_Image(:,:,i,k),sigma_list(j));
            method_grad = 'prewitt';
            [Gx, Gy] = imgradientxy(Smoothed_data,method_grad);
            [Gxx, Gxy] = imgradientxy(Gx);
            [Gyx, Gyy] = imgradientxy(Gy);

            Determinant = Gxx.*Gyy - 2*Gxy;
            Sum_eigenvalues = (Gxx.^2 + Gyy.^2 - 2*Gxy.^2);
           

            List_matrix(:,:,i,k,j) = ((gather(Sum_eigenvalues)));
        end
    end
end

%Extracting the highest determinant obtained 

Cleaned_signal = max(List_matrix,[],5);

%Performing discretization using Triangle method

disp("Performing automated thresholding...")

Discretized_signal = zeros(size(Processed_Image));

for k = 1:l

    Cleaned_signal_temp = real(Cleaned_signal(:,:,:,k));
    
    X = Cleaned_signal_temp(Cleaned_signal_temp>0);
    %Sub-sampling a fraction of the pixels to reduce computation time :
    X = datasample(X,10000);
    Fitted_distribution = fitdist(X,'Gamma');    
    Threshold = Fitted_distribution.icdf(1-Quantile_parameter);
    Discretized_signal(:,:,:,k) = Cleaned_signal_temp > Threshold;
    
    
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
    Position_spots_filtered = Position_spots(Size_spots>1,:);
    List_spots{k}=Position_spots_filtered;
        
end


end

