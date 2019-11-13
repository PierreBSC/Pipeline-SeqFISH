%%%Define a new experiment 
Design_matrix = ["DAPI" "RNA" "RNA" "RNA"];

N_chanel =size(Design_matrix,2);
N_round = size(Design_matrix,1);
N_position = 1;

Design_matrix = repmat(Design_matrix,[1 1 N_position]);


Result_array = cell([size(Design_matrix) N_position]); %%Cell array containing the results of the analysis of each chanel (except DAPI)
Hull_array = cell([size(Design_matrix) N_position]); %%Cell array containing the results of the segmentation based on each chanel (including DAPI)

SNR_table = cell([size(Design_matrix) N_position]);

main_directory = "/home/pbost/Desktop/Image_database/Fat_2/";

temporary_directory = "/Users/Pierre/Desktop/seqFISH_project/Temporary_file/";



%%%RNA spot detection 
for R=1:N_round
    
    if sum(sum(Design_matrix=="RNA"))==0
        break
    end

    Round_directory = strcat(main_directory,"Round_",string(R),"/");
    Round_directory = char(Round_directory);

    for P=1:N_position
    
    Position_directory = strcat(Round_directory,"Position_",string(P),"/");
    Position_directory = char(Position_directory);
            
    RNA_chanel = find(Design_matrix(R,:,1)=="RNA");
    %Loading only RNA data
    RNA_data=LoadImage(Position_directory,true,RNA_chanel); 
    
    RNA_data = Pre_processing(RNA_data); 
    Unfiltered_spots=H_dome_filter(RNA_data);  %%Identification of all possible spots 
    Intensity_threshold = Spot_parameter_optimisation(RNA_data,Unfiltered_spots);
    [Filtered_spots] =Spot_filtering(Unfiltered_spots,Spot_intensity,Intensity_threshold); %%Filtering of the spots
    
    Result_array(R,RNA_chanel,P) = Filtered_spots;
    
    
    %Density based clustering of the spots and export them
    for j=1:size(RNA_data,4)
        segmentation = density_based_segmentation(Filtered_spots{j}(:,1:2))
        Hull_array{R,RNA_chanel(j),P} = segmentation;
    end
    
    %clear RNA_data
    
    end
end






%%%DAPI nuclei segmentation

for R=1:N_round
    
    if sum(Design_matrix(R,:,1)=="DAPI")>1
        disp("Error : several DAPI channels identified") %%If multiple DAPI channel -> error....
        break
    end
    
    Round_directory = strcat(main_directory,"Round_",string(R),"/");
    Round_directory = char(Round_directory);
    
    for P=1:N_position
        Position_directory = strcat(Round_directory,"Position_",string(P),"/");
        Position_directory = char(Position_directory);

        DAPI_chanel = find(Design_matrix(R,:,1)=="DAPI");
        DAPI_data=LoadImage(Position_directory,true,DAPI_chanel); 
        DAPI_data = Pre_processing(DAPI_data,true,true,true); 
        [Nuclei_hull,Morpho_analysis] = Nuclei_identification(DAPI_data,800);
        
        Hull_array{R,DAPI_chanel,P} = Nuclei_hull;

    end
    
        
end



%%%IF analysis 

for R=1:N_round
    
    if sum(sum(Design_matrix=="IF"))==0
        break
    end
    
    Round_directory = strcat(main_directory,"Round_",string(R),"/");
    Round_directory = char(Round_directory);
    
    
    for P=1:N_position
    
    Position_directory = strcat(Round_directory,"Position_",string(P),"/");
    Position_directory = char(Position_directory);
            
    IF_chanel = find(Design_matrix(R,:)=="IF"); %%Select IF images only
    IF_data=LoadImage(Position_directory,true,IF_chanel); %%Load IF images only
    
    IF_data = Pre_processing(IF_data,true,true,false); %Cleaning/scaling of the images
    [IF_objects IF_hull] = IF_analysis(IF_data);  %%Identification of all IF positive cells 
    
    Result_array(P,IF_chanel) = IF_objects;
    
        for j = 1:size(IF_chanel,1)
            Hull_array{R,IF_chanel(j),P} = IF_hull{j};
        end
    end
end




%%%%Alignment of the pictures


if N_round > 1 
        
for R=1:(N_round-1)
    
    if sum(Design_matrix(R,:)=="DAPI")>1
        disp("Error : several DAPI channels identified") %%If multiple DAPI channel -> error....
        break
    end
    
    Ref_round_directory = strcat(main_directory,"Round_",string(R),"/");
    Ref_round_directory = char(Ref_round_directory);
    
    Moving_round_directory = strcat(main_directory,"Round_",string(R+1),"/");
    Moving_round_directory = char(Moving_round_directory);

    
    for P=1:N_position
        
        
        Ref_position_directory = strcat(Ref_round_directory,"Position_",string(P),"/");
        Ref_position_directory = char(Ref_position_directory);
                
        Moving_position_directory = strcat(Moving_round_directory,"Position_",string(P),"/");
        Moving_position_directory = char(Moving_position_directory);

        DAPI_chanel = find(Design_matrix(R,:)=="DAPI");
        
        %%Loading the DAPI staining for the round R and R+1

        Ref_DAPI_data=LoadImage(Ref_position_directory,true,DAPI_chanel);
        Moving_DAPI_data=LoadImage(Moving_position_directory,true,DAPI_chanel);
        
        %%The computation of alignemnt is performed stack-wise : mean value
        %%computed across stacks

        Ref_DAPI_data = mean(Ref_DAPI_data,3);
        Moving_DAPI_data =mean(Moving_DAPI_data,3);
        
        %%Computation of the alignment vector 

        Ref_DAPI_data = Pre_processing(Ref_DAPI_data,true,true,true); 
        Moving_DAPI_data = Pre_processing(Moving_DAPI_data,true,true,true); 
        
        Ref_DAPI_data = imhistmatchn(Ref_DAPI_data,Moving_DAPI_data);
        
        
        transf_vector = Image_alignment(Moving_DAPI_data,Ref_DAPI_data);        
        
        Result_array{R,DAPI_chanel,P} = transf_vector;


    end
    
        
end
end


%%Applying the transformation on the points and on the hulls for round
%%alignment 

Result_array_aligned = cell(size(Result_array));
Hull_array_aligned = cell(size(Hull_array));


Result_array_aligned(1,:,:) = Result_array(1,:,:);
Hull_array_aligned(1,:,:) = Hull_array(1,:,:);


if N_round > 1  
    
    for R=2:N_round
        
    DAPI_chanel = find(Design_matrix(R,:)=="DAPI");
    RNA_chanel = find(Design_matrix(R,:)=="RNA");

        for P=1:N_position
            
         transf_vector = Result_array{R-1,DAPI_chanel,P};
            
            
           for C=1:N_chanel

                %%%To begin with : translation of the convex hull

                temp_hull = Hull_array{R,C,P};
                
                 %%%If no analysis performed on this channel -> no
                 %%%transformation 

                if ~isempty(temp_hull)
                     %%%We apply the transformation to each component of
                     %%%the list
                	temp_hull = cellfun(@transf_vector.transformPointsForward,temp_hull,"UniformOutput",false);
                    Hull_array_aligned{R,C,P} = temp_hull;
                end
                
                %%%Second step : if RNA channel ->
                if Design_matrix(R,C,P)=="RNA"
                    temp_spot_position = Result_array{R,C,P};
                    temp_spot_position(:,1:2) = transf_vector.transformPointsForward(temp_spot_position(:,1:2)) ;
                    Result_array_aligned{R,C,P}  = temp_spot_position;
                end
           end    
        end  
    end
    
end




Selected_channel = 1;
figure, imshow(imadjust(DAPI_data(:,:,1)));

figure, imshow((DAPI_data(:,:,15,1)));
hold on
u = Hull_array{1,1};
for i=1:size(u)
    v = polyshape(u{i});
    plot(v)
end

u = Hull_array{1,2};
for i=1:size(u)
    v = polyshape(u{i});
    plot(v)
end
u = Hull_array{1,4};
for i=1:size(u);
    v = polyshape(u{i});
    plot(v)
end
hold off


figure, imshow((RNA_data(:,:,5,2)));
%colormap("jet")
hold on
    plot(Result_array{1,2}(:,1),Result_array{1,2}(:,2),"*r",'MarkerSize',3)
    plot(Result_array{1,2}(:,2),Result_array{1,2}(:,1),"*y",'MarkerSize',3)
   plot(Result_array{1,5}(:,1),Result_array{1,5}(:,2),"*b",'MarkerSize',2)
hold on
u = Hull_array{1,1};
for i=1:size(u);
    v = polyshape(u{i});
    plot(v,'FaceColor','white','FaceAlpha',0.5)
end

hold off

figure, imshow((RNA_data(:,:,1,2)));
%colormap("jet")
hold on
    plot(Result_array{1,2}(:,2),Result_array{1,2}(:,1),"*r",'MarkerSize',3)

    plot(Result_array{1,1}(:,2),Result_array{1,1}(:,1),"*y",'MarkerSize',3)
    plot(Result_array{1,5}(:,2),Result_array{1,5}(:,1),"*g",'MarkerSize',3)

figure, imshow((mean(IF_data,3)));
hold on
u = Hull_array{1,4};
for i=1:size(u);
    v = polyshape(u{i}(:,[2 1]));
    plot(v,'FaceColor','blue','FaceAlpha',0.5)
end

figure, imshow(imadjust(mean(DAPI_data,3)));
hold on
u = Hull_array{1,3};
for i=1:size(u);
    v = polyshape(u{i}(:,[1 2]));
    plot(v,'FaceColor','blue','FaceAlpha',0.5)
end


    plot(Result_array{1,4}(:,1),Result_array{1,4}(:,2),"*r",'MarkerSize',0.5)

hold off

final_table = spot_assignment(Filtered_spots,Hull_array{1,1});

writetable(table(Result_array{1,4}(:,1:2)),"/Users/Pierre/Desktop/table_2.txt","delimiter","\t");

