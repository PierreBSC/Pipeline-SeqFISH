function [spot_list] = H_dome_filter(Processed_Image,use_GPU,sigma_small,sigma_max,h_meta_parameter,S,N_points_sample)
%Spot detection performed using H-dome filtering
%Implementation inspired from "Quantitative Comparison of Spot Detection
%Methods in Fluorescence Microscopy", IEEE transavtions on medical imaging


if nargin < 2
    use_GPU = true;
    sigma_small = 1; %% Size of the LOG filter, in pixel 
    sigma_max = 5; %% Maximum size of the object in pixel 
    h_meta_parameter = 1; %Don't touch it too much...
    S = 10; %Power of the H-dome matrix, has to be bigger than 3 at least..
    N_points_sample = 20000; %% ~ Number of expected spots * 30
    disp('All needed parameters have not been specified. Running using default parameters')
end



%%If needed sending the Image data to the GPU memory
if use_GPU
    if gpuDeviceCount > 0
        Processed_Image = gpuArray(Processed_Image);
        disp('Using GPU computing for pre-processing...');
    end
end



%%Applying Laplacian of Gaussian filter
kdims = [15 15];
h = fspecial('log',kdims, sigma_small);  %define the Laplacian of Gaussian filter
nImage = size(Processed_Image,1);
mImage = size(Processed_Image,2);
numstack = size(Processed_Image,3);
l = size(Processed_Image,4);

filteredImg = zeros(nImage,mImage,numstack,l,'single');

if use_GPU
    if gpuDeviceCount > 0
        filteredImg = gpuArray(filteredImg);
    end
end



for i=1:l
        filteredImg(:,:,:,i) = imfilter(Processed_Image(:,:,:,i), h);%%Apply LOG filter
end

disp('LOG filter applied')

%clear Processed_Image 

%%Rescaling the image and 
filteredImg = - filteredImg;
filteredImg(filteredImg<0) = 0;

for k = 1:l
    temp_data = gather(filteredImg(:,:,:,k));
    filteredImg(:,:,:,k) = imadjustn(temp_data,stretchlim(temp_data(:)));
end

%%Computing the h-dome value for each channel
h_parameter = zeros(l,1);

for k = 1:l
    temp_data = filteredImg(:,:,:,k);    
    noise = std(gather(temp_data(:)));
    h_parameter(k,1) = (1-noise)*h_meta_parameter;    
end

%%H-dome transformation

H_dome = zeros(nImage,mImage,numstack,l,'single');

if use_GPU
    if gpuDeviceCount > 0
        H_dome = gpuArray(H_dome);
    end
end


for k=1:l
    for z =1:numstack   
        H_dome(:,:,z,k) =  imreconstruct(filteredImg(:,:,z,k)-h_parameter(k),filteredImg(:,:,z,k),4);
        H_dome(:,:,z,k) = filteredImg(:,:,z,k) - H_dome(:,:,z,k);
    end
end 
disp('Morphological Grayscale Reconstruction done')
clear filteredImg;

%%To enhance signal : we put the matrix to the power S

disp('Computation of matrix power S')
H_dome_powered = H_dome.^S; 

disp('Starting Monte-Carlo random sampling')

Sampling_result = cell(1,l);

%%We now sample the pixels using the H_dome_powered probability 
for k=1:l   
    %%We normalise the matrix to sum up to 1
    A = H_dome_powered(:,:,:,k); 
    s = sum(sum(sum(A)));
    A = A./s;   
    A = gather(A);
    
    %%We transform A into a vector and remember the position of each
    %%element
    n_stack = size(H_dome_powered,3);
    
    row_vals = [1:size(A,1)]'*ones(1,size(A,2));  %all x values
    row_vals = row_vals(:);
    row_vals = repmat(row_vals,n_stack,1);
    
    col_vals = ones(size(A,2),1)*[1:size(A,2)];  %all y values
    col_vals = col_vals(:);
    col_vals = repmat(col_vals,n_stack,1);

    stack_vals = repelem(1:n_stack,size(A,1)*size(A,2)).';
    
    A = A(:);
    
    %%We now sample :
    sampling = randsample(1:size(A),N_points_sample,true,A);
    Position = [row_vals(sampling),col_vals(sampling),stack_vals(sampling)];

    Sampling_result{k} = Position;
    disp('....')
end

disp('Sampling done')

spot_list = cell(l,1);

disp('Starting Mean shift clustering')
%%We now cluster the sampled pixels using Mean Shift algorithm


for k=1:l
    
    %%Clustering by itself
    [Final_clusters,Final_clustering] = Mean_shift(Sampling_result{k},sigma_small*3,20);
    size(Final_clusters)
    
    count_cluster = tabulate(Final_clustering);
        
    grouped_spots = splitapply( @(x){x}, Sampling_result{k}, (Final_clustering));
    
    %%Computing the determinant of the spatial variance for each cluster
    spot_determinant = cellfun(@(x) {det(cov(x(:,1:2)))},grouped_spots);
    spot_determinant = cell2mat(spot_determinant);
    
    %Computing the maximal determinant value 
    det_threshold = ((sigma_max^2)^2);
    
     %Filtering the spots based on the different values  
    Selected_spots = spot_determinant<det_threshold & count_cluster(:,2)>5;

    
    spot_list{k} = Final_clusters(Selected_spots,[2 1 3]);

end

end

