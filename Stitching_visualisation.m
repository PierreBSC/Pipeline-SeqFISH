function [] = Stitching_visualisation(Analysis_result,Parameters,R,Channel)
%Simple function to visualy check the quality of the image tiling

if nargin < 3
    R = 1;
    Channel = 1;
end


if ~isfield(Analysis_result,'Global_stitching')
    disp("Perform image stitching step before running this function")
    return
end

M = -Analysis_result.Global_stitching;
X_min = min(M(:,1));
Y_min = min(M(:,2));
X_max = max(M(:,1));
Y_max = max(M(:,2));

%Adding extra values to capture the whole pocute 

X_min = X_min - 1000;
Y_min = Y_min - 1000;
X_max = X_max + 1000;
Y_max = Y_max + 1000;


%Creating the (empty) merged image)
Mixed_image = zeros(X_max-X_min,Y_max-Y_min,Parameters.N_position);

%Creating the spatial reference
panorama_view = imref2d([X_max-X_min,Y_max- Y_min],[X_min X_max], [Y_min Y_max]);


for k = (1:Parameters.N_position)
    
    T = affine2d([1 0 0 ; 0 1 0 ; -Analysis_result.Global_stitching(k,1) ,-Analysis_result.Global_stitching(k,2) 1]);
    
    Round_directory = strcat(Parameters.Image_directory,"/Round_",string(R),"/");
	Round_directory = char(Round_directory);
            
	Position_directory = strcat(Round_directory,"/Position_",string(k),"/");
	Position_directory = char(Position_directory);

    
    Image = LoadImage(Position_directory,true,Channel); 
    Image = Pre_processing(Image,Parameters.use_GPU,Parameters.Substack,Parameters.Stack_min,Parameters.Stack_max,Parameters.perform_background_removal,Parameters.background_sigma_parameter,false,0); 
	Image = mean(Image,3);
    warped_Image = imwarp(Image,T,'OutputView',panorama_view);

    Mixed_image(:,:,k) = warped_Image;
    
end

figure, imshow(imadjust(max(Mixed_image,[],[3])),[])




end

