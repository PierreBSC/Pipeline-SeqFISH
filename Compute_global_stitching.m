function [Analysis_result] = Compute_global_stitching(Analysis_result,Parameters,R,Channel)
%General function to compute the stitching vector for all positions 
%This function uses a phase correlation based method to efficiently compute
%each individual vector and then use a graph based approach to compute the
%final translation vector 

if nargin < 4
    Channel = 1;
    R = 1;
end 


if Parameters.use_GPU
    if gpuDeviceCount > 0
        d = gpuDevice;
        fprintf('\nWill use GPU computing for Image stitching steps: %s.\n',d.Name)
    end
end



%%Creating the field of the Analysis_result object : list of vectors to
%%apply to each

Analysis_result.Global_stitching = zeros(Parameters.N_position,2);

%Creating the three matrix used : 
% M_x and M_y elements (i,j) correponds to the x and y component of the
% translation vector to align image i to image j
% M_SNR elements (i,j) : quality of the peaks 

M_x = zeros(Parameters.N_position,'single');
M_y = zeros(Parameters.N_position,'single');
M_SNR = zeros(Parameters.N_position);


%First step : computing each possible translation vector for each pair of
%positions 

for i = 1:Parameters.N_position
    parfor j =1:Parameters.N_position
        if i~=j
            
            Round_directory = strcat(Parameters.Image_directory,"/Round_",string(R),"/");
            Round_directory = char(Round_directory);
            
            Position_directory_i = strcat(Round_directory,"/Position_",string(i),"/");
            Position_directory_i = char(Position_directory_i);
            
            Position_directory_j = strcat(Round_directory,"/Position_",string(j),"/");
            Position_directory_j = char(Position_directory_j);

            Ref_Image = LoadImage(Position_directory_i,true,Channel); 
            Ref_Image = Pre_processing(Ref_Image,Parameters.use_GPU,Parameters.Substack,Parameters.Stack_min,Parameters.Stack_max,Parameters.perform_background_removal,Parameters.background_sigma_parameter,true,0.001); 
            Ref_Image = mean(Ref_Image,3);
            Ref_Image = imadjust(imgaussfilt(Ref_Image,2));

            
            Moving_Image = LoadImage(Position_directory_j,true,Channel); 
            Moving_Image = Pre_processing(Moving_Image,Parameters.use_GPU,Parameters.Substack,Parameters.Stack_min,Parameters.Stack_max,Parameters.perform_background_removal,Parameters.background_sigma_parameter,true,0.001); 
            Moving_Image = mean(Moving_Image,3);
            Moving_Image = imadjust(imgaussfilt(Moving_Image,2));
            
            
            if Parameters.use_GPU
                if gpuDeviceCount > 0
                    Moving_Image = gpuArray(Moving_Image);
                    Ref_Image = gpuArray(Ref_Image);
                end
            end

            [Stitching_transformation,Intensity_peak] = get_stitching_vector(Moving_Image,Ref_Image,"Inversed");

            M_SNR(i,j) = Intensity_peak;
            Stitching_transformation = Stitching_transformation.T;
            M_x(i,j) = gather(Stitching_transformation(3,1));
            M_y(i,j) = gather(Stitching_transformation(3,2));
            
        end
    end
end

%%Performing automated thresholding to identify the real translocation

List_SNR = triu(M_SNR);
List_SNR = List_SNR(:);
List_SNR = List_SNR(List_SNR~=0);
List_SNR = sort(List_SNR,'descend');

%Estimating the upper and lower bound of the threshold

%At least P-1 connections between the images to get a connected graph
Threshold_max = List_SNR(Parameters.N_position - 1);

%If images are coming from a square grid we can estimate the number of
%significant vectors
Threshold_min = List_SNR(round(2*Parameters.N_position-2*sqrt(Parameters.N_position)));

%Performing the thresholding by itself
[count binning_values] = hist(List_SNR,20);
Threshold = otsuthresh(count);
Threshold = Threshold*max(binning_values);

%%Making some checks on the threshold value :
if Threshold > Threshold_max
    Threshold = Threshold_max
end

if Threshold < Threshold_min
    Threshold = Threshold_min
end


%Creating the graph
M_SNR_thresholded = M_SNR >Threshold ;
Graph_SNR = graph((M_SNR_thresholded+transpose(M_SNR_thresholded))/2)

%Is the threshold value giving a unique connected component ?

[~,comp_count]=conncomp(Graph_SNR);
N_significant_conn_comp = sum(comp_count>1);
New_threshold = [];

In_range = true;
Unique_comp_reached = false;

if N_significant_conn_comp==1
    Unique_comp_reached = true;
end


List_SNR_in_range = List_SNR( List_SNR > Threshold_min & List_SNR < Threshold);
k=1;

while  In_range & ~Unique_comp_reached
    
    Temp_threshold = List_SNR_in_range(k);
    Temp_graph = (M_SNR>Temp_threshold);
    Temp_graph = graph((Temp_graph + transpose(Temp_graph))/2);
    [~,comp_count_temp]=conncomp(Temp_graph)
    N_significant_conn_comp_temp = sum(comp_count_temp>1);
    
    if N_significant_conn_comp_temp == 1
        Unique_comp_reached = true;
    end
    k = k +1;
    
    if k > size(List_SNR_in_range)
       In_range = false; 
    end
    
    Graph = Temp_graph;

end


List_translation = [];

for k = 1:Parameters.N_position

    temp_path = shortestpath(Graph_SNR,k,1);
    d = size(temp_path,2)-1;
    List_dist_X = [];
    List_dist_Y = [];

    for i= 1:d
        List_dist_X = [List_dist_X M_x(temp_path(i),temp_path(i+1)) ];
        List_dist_Y = [List_dist_Y M_y(temp_path(i),temp_path(i+1)) ];
    end
    Mov_X = sum(List_dist_X);
    Mov_Y = sum(List_dist_Y);
    List_translation  = [List_translation ; [Mov_X,Mov_Y]];

end

%Checking the final results
plot(Graph_SNR)
Analysis_result.Global_stitching = -List_translation;

end

