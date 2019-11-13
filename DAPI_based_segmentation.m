function [Analysis_result] = DAPI_based_segmentation(Analysis_result,Parameters)
%Simple function to perform cell segmentation based on nuclear DAPI
%staining
%The script is based on image binarisation followed by Euclidean Distance
%Map and Watershed transformed


%First step : create a new element inside the Analysis_result to store the
%results of segmentation

Analysis_result.Nuclei_segmentation = cell(Parameters.N_round,Parameters.N_channel,Parameters.N_position);


for R = 1:Parameters.N_round
    
    %If no DAPI channel at this round : go to next round...
    if sum(sum(Parameters.Matrix_design{R,:}=="DAPI"))==0
        continue
    end
    
    %If more than one DAPI channel : error !
    if sum(sum(Parameters.Matrix_design{R,:}=="DAPI"))>1
        dips('ERROR : more than one DAPI channel in the design matrix !')
        break
    end
    

    
    %Defining the Round directory
    Round_directory = strcat(Parameters.Image_directory,'/Round_',string(R),'/');
    Round_directory = char(Round_directory);

    for P = 1:Parameters.N_position
    
    disp(strcat('Round ',num2str(R),' Position ',num2str(P))) 
    
    %Defining the Position subdirectory
    Position_directory = strcat(Round_directory,'/Position_',string(P),'/');
    Position_directory = char(Position_directory);
    
    %Which are the RNA channels ?
    DAPI_channel = find(Parameters.Matrix_design{R,:}=="DAPI");
    
    %Loading only RNA data
    DAPI_data=LoadImage(Position_directory,true,DAPI_channel); 
    
    if Parameters.use_GPU
        if gpuDeviceCount > 0 
            DAPI_data = gpuArray(DAPI_data);
            disp('Using GPU computing for pre-processing...');
        end
    end

    
    %Cleaning of the data
    DAPI_data = Pre_processing(DAPI_data,Parameters.use_GPU,Parameters.Substack,Parameters.Stack_min,Parameters.Stack_max,Parameters.perform_background_removal,Parameters.background_sigma_parameter,Parameters.perform_intensity_adjustment,Parameters.tolerance); 
    DAPI_data = max(squeeze(DAPI_data),[],3);
    
    %Ajusting the intensity : 1% tolerance is usually enough
    DAPI_data = imadjust(DAPI_data);
    
    %Binarisation of the pictures and filling the wholes
    %To avoid over-segmentation : use of image closing
    
    DAPI_data_bin = imbinarize(DAPI_data); %%Otsu's thresholding
    DAPI_data_bin = imclose(DAPI_data_bin,strel('disk',4));
    DAPI_data_bin = imfill(DAPI_data_bin,'holes'); %%Otsu's thresholding

    %%Segmentation by itself
    
    D = -bwdist(~DAPI_data_bin);  %%Distance map
    D(~DAPI_data_bin) = Inf;
 
    mask = imextendedmin(D,0); 
    D2 = imimposemin(D,mask); %%Identifying h-bassin
 
    Watershed_analysis = watershed(D2);
    Watershed_analysis(~DAPI_data_bin) = 0;

Nuclei_hull_list = regionprops(Watershed_analysis,"ConvexHull");
Nuclei_hull_list = struct2cell(Nuclei_hull_list).';

Nuclei_morphological_properties = regionprops(Watershed_analysis,"Eccentricity","Orientation","Solidity");
Nuclei_morphological_properties = struct2table(Nuclei_morphological_properties);

Nuclei_size = struct2table(regionprops(Watershed_analysis,"Area"));

Nuclei_hull_list_filtered = Nuclei_hull_list(Nuclei_size.Variables>=min_size);
Nuclei_morphological_properties_filtered = Nuclei_morphological_properties(Nuclei_size.Variables>=min_size,:);
    end
end


end

