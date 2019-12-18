function [Cleaned_segmentation] = Remove_segmentation_overlap(Segmentation_results,X_size,Y_size)
%Clean the segmentation results to remove overlaping area 
% This subroutine uses the results of a segmentaiton under the shape of a
%pixel table

N_cells = max(Segmentation_results(:,3));

%First step : checking that we do have overlap 
Segmentation_results(Segmentation_results==0)=1;

Count_pixels = tabulate((sub2ind([X_size Y_size],Segmentation_results(:,1),Segmentation_results(:,2))));
Count_pixels = tabulate(Count_pixels(:,2)>1);

%If there is indeed overlap then
if size(Count_pixels,1)>1
    
    ED_map_list = cell(N_cells,1);
    
    %%Constructing ED map
    for k = 1:N_cells
        list_pixels_temps = Segmentation_results(Segmentation_results(:,3)==k,1:2);
        M = zeros(X_size,Y_size);
        M(sub2ind([X_size Y_size],list_pixels_temps(:,1),list_pixels_temps(:,2))) = 1;
        M = bwdist(~M);
        M = M./sum(M(:));
        ED_map_list{k} = M;
    end
    
    Polyshape_list = cell(N_cells,1);
    %Creating polyshape objects for each cell
    for k = 1:N_cells 
        X = Segmentation_results(Segmentation_results(:,3)==k,1:2);
        Final_boundary = X(boundary(X(:,1:2)),:);
        Polyshape_list{k} = polyshape(Final_boundary);
    end
    
    %Which cells are overlaping ?
    
    Overlap_matrix = zeros(N_cells);
    
    for i=1:N_cells
        for j=1:N_cells
            s = overlaps(Polyshape_list{i},Polyshape_list{j});
            Overlap_matrix(i,j) = s;
        end
    end
    
    Overlap_matrix = triu(Overlap_matrix,1);
    [x ,y] = ind2sub([N_cells, N_cells],find(Overlap_matrix));
    List_interaction = [x,y];
    N_interactions = size(List_interaction,1);
    
    %Now we can correct the overlaps in an iterative manner
    for l=1:N_interactions
        
        i = List_interaction(l,1);
        j = List_interaction(l,2);
        
        if overlaps(Polyshape_list{i},Polyshape_list{j})
            %Attribution of the pixels to a given cell by looking at the
            %arg.max value of the normalised Euclidean distance map
            M=sign(ED_map_list{i}-ED_map_list{j});
            %Changing the two polyshape objects 
            [x_i,y_i] = ind2sub([X_size Y_size],find(M==1));
            X_i = [x_i,y_i];
            [x_j,y_j] = ind2sub([X_size Y_size],find(M==(-1)));
            X_j = [x_j,y_j];

            new_boundary_i = X_i(boundary(X_i),:);
            new_boundary_j = X_j(boundary(X_j),:);
            Polyshape_list{i} = polyshape(new_boundary_i);
            Polyshape_list{j} = polyshape(new_boundary_j);
            
        end
        
    end

     %Now we can export the results of the experiment
    Cleaned_segmentation = [];
    for k=1:N_cells
        X = Polyshape_list{k};
        M = poly2mask(X.Vertices(:,1),X.Vertices(:,2),Y_size,X_size);
        [x_temp, y_temp] = ind2sub([X_size Y_size],find(M==1));
        temp_results = [y_temp, x_temp, repmat(k,size(x_temp,1),1)];
        Cleaned_segmentation = [Cleaned_segmentation;temp_results];
    end
    
    
    
end

if size(Count_pixels,1)==1
    Cleaned_segmentation = Segmentation_results
end


end

