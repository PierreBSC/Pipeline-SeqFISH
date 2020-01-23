function Parameters = Define_data_path(Parameters)
% Ask user to define path containing the data

if isempty(Parameters) 
    path = '';
else
    path = Parameters.Image_directory;    
end
    
Image_directory = uigetdir(path,'Provide the path to the Image Folder');
if Image_directory ~= 0 
    Parameters.Image_directory = Image_directory;
end


