function [Analysis_result,Parameters] = Create_experiment()
%First function of the pipeline 
%Create the different variables, load the settings for the analysis

%Geting the main image directory
Image_directory = uigetdir('title','Provide the path to the Image Folder');


%Geting the parameter through a user Interface
Parameters = Define_parameters();


%Geting the matrix design and loading it

[x , y] = uigetfile('*','Provide the path to the Matrix design file');
Matrix_design_path = strcat(y,x);
Matrix_design = readtable(Matrix_design_path,'delimiter','\t');

Parameters.Matrix_design = Matrix_design;
Parameters.N_round = size(Matrix_design,1);
Parameters.N_channel = size(Matrix_design,2);
Parameters.Image_directory = Image_directory;

%Creating the variable containing all the results of analysis

Analysis_result = struct;

end

