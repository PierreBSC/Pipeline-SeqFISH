function [Analysis_result,Parameters] = Create_experiment()
Analysis_result = {};
Parameters = {};

%First function of the pipeline 
%Create the different variables, load the settings for the analysis

% Aks user for the main image directory
Image_directory = uigetdir('title','Provide the path to the Image Folder');
if Image_directory == 0; return; end
Parameters.Image_directory = Image_directory;

%Geting the matrix design and loading it
Parameters = Define_data_organization(Parameters);
if isempty(Parameters); return; end


% Geting the parameter for spot detection through a user Interface
Parameters = Define_parameters_spot_detection(Parameters);
if isempty(Parameters); return; end

% Adding parameters for cell segmentation
Parameters = Define_parameters_cell_segmentation(Parameters);
if isempty(Parameters); return; end

%Creating the variable containing all the results of analysis
Analysis_result = struct;

end

