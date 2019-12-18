function [Analysis_result,Parameters] = Spot_based_segmentation(Analysis_result,Parameters)
%Generic function that performs cell segmentation using spots
%Requires an Analysis_result and Parameters objects to already exist

%Checking that the spots detection stage has already been performed

if ~isfield(Analysis_result,'Spot_analysis_raw')
    disp("Perform spot detection step before running this function")
    return
end

if Parameters.Perform_on_filtered_spots
    if ~isfield(Analysis_result,'Spot_analysis_filtered')
        Parameters.Perform_on_filtered_spots =false;
    end
end

%Now we can perform cell segmentation by itself
%Creating a new field in the Analysis_result object

Analysis_result.Spot_based_segmentation = cell(Parameters.N_position,1);

%Extracting spot array

if Parameters.Perform_on_filtered_spots
    Spot_array = Analysis_result.Spot_analysis_filtered;
end

if ~Parameters.Perform_on_filtered_spots
    Spot_array = Analysis_result.Spot_analysis_raw;
end

%%Asking which channel to use for segmentation
Segmentation_design = Select_segmentation_genes(Parameters);

%%We need to know what is the size of the pictures : just load one
%%picture...
Round_directory = strcat(Parameters.Image_directory,"/Round_1/");
Round_directory = char(Round_directory);
Position_directory = strcat(Round_directory,"/Position_1/");
Position_directory = char(Position_directory);
Example_data=LoadImage(Position_directory,true,1); 
X_size = size(Example_data,1);
Y_size = size(Example_data,2);


%Performing the Diffusion based clustering method on each RNA channel on
%each Position and for each Round

Spot_clusters = cell(Parameters.N_round,Parameters.N_channel,Parameters.N_position);


disp("Perfoming diffusion based clustering and clustering filtering")
for R = 1:Parameters.N_round
    
    %If no RNA channel at this round : go to next round...
    if sum(sum(Parameters.Matrix_design{R,:}=="RNA"))==0
        continue
    end
    
    %What are the RNA channels for this round ?
    Segmentation_channels = find(Parameters.Matrix_design{R,:}=="RNA" & Segmentation_design{R,:}=="Use");

    for P = 1:Parameters.N_position
    
    disp(strcat("Round ",num2str(R)," Position ",num2str(P))) 
    
        for k=1:size(Segmentation_channels,2)
            
            temp_spots = Spot_array{R,Segmentation_channels(1,k),P};
            %%Key step :Computing NJW clustering
            temp_clustering = Diffusion_map(temp_spots,Parameters.N_neighbors,Parameters.N_comp,Parameters.T);
            
            
            %Cleaning of the different spots  clusters 
            %Beginning with the number of spots : 10 at least is a good
            %start
            
            distribution_cluster = tabulate(temp_clustering);
            list_cluster =unique(temp_clustering);
                        
            %%Second filter : are the c
            %What is the spot density in each cell/cluster ?
            
            Lambda = size(temp_spots,1)/((max(temp_spots(:,2))-min(temp_spots(:,2)))*(max(temp_spots(:,1))-min(temp_spots(:,1)))) ;
            Volume_cluster = [];
            P_value_poisson =[];
            
            for i=1:size(list_cluster,1)
                
                %Checking for the density of the cluster using a simple
                %Homogenous Poisson Point Process model 
                if sum(temp_clustering==i)>2
                    [~,volume] = convhull(temp_spots(temp_clustering==i,1:2));
                    Poisson_distribution = makedist('Poisson','lambda',Lambda*volume); 
                    p = 1-cdf(Poisson_distribution,distribution_cluster(i,2));
                    P_value_poisson = [P_value_poisson p];
                end
                
                if sum(temp_clustering==i)<=2
                    volume = 0;
                    P_value_poisson = [P_value_poisson 1];
                end
                


            end
            
            
            %%Good quality spot clusters : more than Minimal_spot_cells spots and high
            %%spot density
            %Minimal_spot_cells value should be around 10 for high quality
            %cells....
            selected_clusters = distribution_cluster(distribution_cluster(:,2)>Parameters.Minimal_spot_cells & P_value_poisson.' < 0.0001,1);
            temp_spots_filtered = temp_spots(ismember(temp_clustering,selected_clusters),:);
            temp_clustering_filtered = temp_clustering(ismember(temp_clustering,selected_clusters));
            Final_table = [temp_spots_filtered temp_clustering_filtered];
            
            
           Spot_clusters{R,Segmentation_channels(1,k),P} = Final_table;       
           
        end
    end
    
end

%%%Now we can aggregate the clusters across genes (i.e across Rounds and
%%%Channels)
%%To do so : compute the estimated shape through density estimation of the
%%spots : way more robust than convex hull estimation

disp("Estimating cell shape through density computation")

Cluster_density_map = cell(size(Spot_clusters));


for R = 1:Parameters.N_round
    
    %If no RNA channel at this round : go to next round...
    if sum(sum(Parameters.Matrix_design{R,:}=="RNA"))==0
        continue
    end
    
    %What are the RNA channels for this round ?
    Segmentation_channels = find(Parameters.Matrix_design{R,:}=="RNA");

    for P = 1:Parameters.N_position
    
    disp(strcat("Round ",num2str(R)," Position ",num2str(P))) 
    
        for k=1:size(Segmentation_channels,2)
            
            temp_spots = Spot_clusters{R,Segmentation_channels(1,k),P};
            
            if size(temp_spots,1)>0
             %Re-assigning values of the clusters
            [~, ~, temp_spots(:,4)] = unique(temp_spots(:,4));
            temp_spots(:,4) = discretize(temp_spots(:,4),size(unique(temp_spots(:,4)),1));
            temp_density = splitapply(@(x) {histcounts2(x(:,2),x(:,1),linspace(1,X_size,50),linspace(1,Y_size,50))./size(x,1)},temp_spots(:,1:2),temp_spots(:,4));
            
            Spot_clusters{R,Segmentation_channels(1,k),P} = temp_spots;
            Cluster_density_map{R,Segmentation_channels(1,k),P} = temp_density;

            end
            
        end
    end
    
end

%%Last step : for each Position : aggregating the enveloppes 

disp("Aggregating the spot clusters")


for P = 1:Parameters.N_position
    
    temp_density_list = Cluster_density_map(:,:,P);
    temp_density_list_aggregated = [];
    Cluster_ID = [];
    
    %Collecting all enveloppes from all rounds/channels
    
    for R = 1:Parameters.N_round
        
            for k=1:Parameters.N_channel
                
                temp_density_list_aggregated = [temp_density_list_aggregated ; temp_density_list{R,k}];
                 %Where does each cluster come from and wg
                Cluster_ID = [Cluster_ID [repmat(R,1,size(temp_density_list{R,k},1));repmat(k,1,size(temp_density_list{R,k},1)); (1:size(temp_density_list{R,k},1))] ];
            end
    end
    
    %%Creating a similarity matrix between the enveloppes
    Similarity_matrix = zeros(size(temp_density_list_aggregated,1),size(temp_density_list_aggregated,1));
    
    for i=1:size(Similarity_matrix,1)
        
        for j=1:size(Similarity_matrix,2)
            
            %Computing the overlap score between the two distribution 

            Overlap_score = sqrt(temp_density_list_aggregated{i}.*temp_density_list_aggregated{j});
            Similarity_matrix(i,j) = sum(Overlap_score(:));
            
        end
    end
    
    Aggregation_graph = graph(Similarity_matrix > Parameters.Graph_Threshold_overlap); %
	Aggregation_graph = rmedge(Aggregation_graph, 1:numnodes(Aggregation_graph), 1:numnodes(Aggregation_graph)); %Removing self loops
    Meta_clustering = conncomp(Aggregation_graph); %Detecting connected components 
    N_meta_clusters = size(unique(Meta_clustering),2);
    
    %%Computing the final enveloppes :
    Meta_clusters_points = [];

    for k = 1:N_meta_clusters
        
        Cluster_ID_selected = Cluster_ID(:,Meta_clustering==k);
        temp_point_aggregate = [];
        
        for j=1:size(Cluster_ID_selected,2)
            X = Spot_clusters{Cluster_ID_selected(1,j),Cluster_ID_selected(2,j),P};
            X = X(X(:,4)==Cluster_ID_selected(3,j),:);
            temp_point_aggregate = [temp_point_aggregate;X(:,1:3)];
        end
        
                
        %Computing boundary box
        [boundary_x,boundary_y] = boundingbox(polyshape(temp_point_aggregate(:,1:2)));
        %Extending it
        boundary_x(1) = boundary_x(1) - 2*Parameters.Bandwidth_parameter;
        boundary_y(1) = boundary_y(1) - 2*Parameters.Bandwidth_parameter;
        boundary_x(2) = boundary_x(2) + 2*Parameters.Bandwidth_parameter;
        boundary_y(2) = boundary_y(2) + 2*Parameters.Bandwidth_parameter;
        %Checking it is not outside of the image boundary
        boundary_x(1) = max(0,boundary_x(1));
        boundary_y(1) = max(0,boundary_y(1));
        boundary_x(2) = min(X_size,boundary_x(2));
        boundary_y(2) = min(Y_size,boundary_y(2));
        
        
        %Computing local density using kernel estimator
        [grid_x, grid_y] = meshgrid(round(boundary_x(1)):round(boundary_x(2)),round(boundary_y(1)):round(boundary_y(2)));
        estimated_density = ksdensity(temp_point_aggregate(:,1:2),[grid_x(:),grid_y(:)],'Bandwidth',Parameters.Bandwidth_parameter);
        estimated_density = reshape(estimated_density,size(grid_y));
        
        %Binarization of the density 
        estimated_density = imadjust(estimated_density);
        estimated_density_bin = imbinarize(estimated_density);
        
        %Extracting the biggest connected component
        estimated_density_bin = imfill(estimated_density_bin,'holes'); %Removing holes
        Conn_components = bwconncomp(estimated_density_bin);
        Conn_components = Conn_components.PixelIdxList;
        Conn_components_size = cellfun(@(x) size(x,1),Conn_components);
        Conn_components = Conn_components{Conn_components_size==max(Conn_components_size)};
        %Getting the final positions 
        Conn_components_points = [grid_x(Conn_components),grid_y(Conn_components)];
        Conn_components_points = [Conn_components_points, repmat(k,size(Conn_components_points,1),1)];
        Meta_clusters_points = [Meta_clusters_points; Conn_components_points];
        
        Final_boundary = Conn_components_points(boundary(Conn_components_points(:,1:2)),:);
    end
        
    Analysis_result.Spot_based_segmentation{P,1} = Meta_clusters_points;
    
end
       
for P=1:Parameters.N_position
        Analysis_result.Spot_based_segmentation{P,1} = Remove_segmentation_overlap(Analysis_result.Spot_based_segmentation{P,1},X_size,Y_size);
end

disp("Spot based segmentation finished !")

end

