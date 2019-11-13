function [Nuclei_hull_list_filtered,Nuclei_morphological_properties_filtered] = Nuclei_identification(DAPI_data,min_size)
%Extraction of the core of the nuclei and conversion into a convex hull
%object with Npoints_hull

if nargin <2
    min_size = 50;
end

Mean_image  = (mean(DAPI_data,3));

Mean_image_th = imbinarize(Mean_image); %%Otsu's thresholding

D = -bwdist(~Mean_image_th);  %%Distance map
D(~Mean_image_th) = Inf;
 
mask = imextendedmin(D,2); 
D2 = imimposemin(D,mask); %%Identifying h-bassin
 
Watershed_analysis = watershed(D2);
Watershed_analysis(~Mean_image_th) = 0;

Nuclei_hull_list = regionprops(Watershed_analysis,"ConvexHull");
Nuclei_hull_list = struct2cell(Nuclei_hull_list).';

Nuclei_morphological_properties = regionprops(Watershed_analysis,"Eccentricity","Orientation","Solidity");
Nuclei_morphological_properties = struct2table(Nuclei_morphological_properties);

Nuclei_size = struct2table(regionprops(Watershed_analysis,"Area"));

Nuclei_hull_list_filtered = Nuclei_hull_list(Nuclei_size.Variables>=min_size);
Nuclei_morphological_properties_filtered = Nuclei_morphological_properties(Nuclei_size.Variables>=min_size,:);
disp("Nuclei analysis done")
end

