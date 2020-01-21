P =1;
R=1;

Round_directory = strcat(Parameters.Image_directory,'/Round_',string(R),'/');
Round_directory = char(Round_directory);

disp(strcat('Round ',num2str(R),' Position ',num2str(P))) 

%Defining the Position subdirectory
Position_directory = strcat(Round_directory,'/Position_',string(P),'/');
Position_directory = char(Position_directory);

%Loading only RNA data
RNA_data=LoadImage(Position_directory,true,RNA_channel); 
    
%Cleaning of the data
RNA_data = Pre_processing(RNA_data,Parameters.use_GPU,Parameters.Substack,Parameters.Stack_min,Parameters.Stack_max,Parameters.perform_background_removal,Parameters.background_sigma_parameter,Parameters.perform_intensity_adjustment,Parameters.tolerance); 
RNA_data = RNA_data(:,:,1,:);

Unfiltered_spots = Multiscale_filter(RNA_data,Parameters.use_GPU,Parameters.sigma_small,Parameters.N_scales,Parameters.T_offset);  

temporary_directory = '~/seqFISH_Temporary_file/';
[Intensity_threshold Spot_intensity] = Spot_parameter_optimisation(RNA_data,Unfiltered_spots,Parameters.use_GPU,temporary_directory,strcat(pwd,'/Scan_test.R'),Parameters.sigma_small); %Computing the value intensity threshold
Filtered_spots = Spot_filtering(Unfiltered_spots,Spot_intensity,Intensity_threshold); %%Filtering of the spots

K = 4 ;

X = Unfiltered_spots{K};
Y = Filtered_spots{K};
imshow(RNA_data(:,:,1,K),[]);
hold on
scatter(X(:,1),X(:,2),'r','filled')
scatter(Y(:,1),Y(:,2),'b','filled')

%%Checking for overlapp across rounds
P=1

R = 1;
Round_directory = strcat(Parameters.Image_directory,'/Round_',string(R),'/');
Round_directory = char(Round_directory);
Position_directory = strcat(Round_directory,'/Position_',string(P),'/');
Position_directory = char(Position_directory);


DAPI_round_1 = LoadImage(Position_directory,true,1); 
DAPI_round_1 = Pre_processing(DAPI_round_1,Parameters.use_GPU,Parameters.Substack,Parameters.Stack_min,Parameters.Stack_max,Parameters.perform_background_removal,Parameters.background_sigma_parameter,Parameters.perform_intensity_adjustment,0.01); 
DAPI_round_1 =DAPI_round_1(:,:,1);  
  
R = 2;
Round_directory = strcat(Parameters.Image_directory,'/Round_',string(R),'/');
Round_directory = char(Round_directory);
Position_directory = strcat(Round_directory,'/Position_',string(P),'/');
Position_directory = char(Position_directory);


DAPI_round_2 = LoadImage(Position_directory,true,1); 
DAPI_round_2 = Pre_processing(DAPI_round_2,Parameters.use_GPU,Parameters.Substack,Parameters.Stack_min,Parameters.Stack_max,Parameters.perform_background_removal,Parameters.background_sigma_parameter,Parameters.perform_intensity_adjustment,0.01); 
DAPI_round_2 =DAPI_round_2(:,:,1);  

R = 3;
Round_directory = strcat(Parameters.Image_directory,'/Round_',string(R),'/');
Round_directory = char(Round_directory);
Position_directory = strcat(Round_directory,'/Position_',string(P),'/');
Position_directory = char(Position_directory);


DAPI_round_3 = LoadImage(Position_directory,true,1); 
DAPI_round_3 = Pre_processing(DAPI_round_3,Parameters.use_GPU,Parameters.Substack,Parameters.Stack_min,Parameters.Stack_max,Parameters.perform_background_removal,Parameters.background_sigma_parameter,Parameters.perform_intensity_adjustment,0.01); 
DAPI_round_3 =DAPI_round_3(:,:,1);  


M_1 = normxcorr2(DAPI_round_1,DAPI_round_2);
imshow(M_1,[]), colormap("jet")
[X_1 Y_1] = find(M_1==max(M_1(:)));
X_1 = X_1 - size(DAPI_round_1,1);
Y_1 = Y_1 - size(DAPI_round_1,2);
amplitude_X_1 = X_1/size(DAPI_round_1,1);
amplitude_Y_1 = Y_1/size(DAPI_round_1,1);




M_2 = normxcorr2(DAPI_round_2,DAPI_round_3);
imshow(M_2,[]), colormap("jet")
[X_2 Y_2] = find(M_2==max(M_2(:)));
X_2 = X_2 - size(DAPI_round_1,1);
Y_2 = Y_2 - size(DAPI_round_1,2);
amplitude_X_2 = X_2/size(DAPI_round_1,1)
amplitude_Y_2 = Y_2/size(DAPI_round_1,1)

