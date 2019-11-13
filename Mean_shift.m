function [Final_clusters,Final_clustering] = Mean_shift_2(X,bandwidth,max_iter)
%Mean shift clustering function 
%Mean shift clustering based on Gaussian kernel computed with Matrix
%%Dscription of the global approach : A review of mean-shif algorithms for
%%clustering

[Seed_values Seed_correspondance Attributed_seed]= unique(X,'row');
N_seeds = size(Seed_values,1);

%%We need to compute a Kd-tree to increase computation speed
%We do it with all the spots and not only the seeds

KD_Tree = KDTreeSearcher(X);

%What is the max distance that will be take into account ? (only true for
%gaussian kernel)
max_distance = bandwidth*sqrt(-log(1e-3));
%What is the minimal distance so that two modes are distinct ?
epsilon_distance = bandwidth/2;

List_modes = cell(size(Seed_values,1),1);

parfor k=1:N_seeds
    
    mean_pos = Seed_values(k,:);
    completed_iterations = 1;
    keep_climb_gradient = true;
    
    while keep_climb_gradient
    
        %Computation fo gaussian kernel
        
        old_mean_pos = mean_pos;  % save the old mean position
        temp_neighours = rangesearch(KD_Tree,old_mean_pos,max_distance);
        temp_neighours = X(temp_neighours{1},:);
                
        Gaussian_weights = exp(-sum((temp_neighours - mean_pos).^2,2) /bandwidth^2);
        mean_pos = sum(Gaussian_weights.*temp_neighours,1)/sum(Gaussian_weights);
    
        %%Has the method converged ?
        if norm(mean_pos-old_mean_pos) < epsilon_distance | completed_iterations == max_iter
             keep_climb_gradient = false ;    
        end
        
       %What happens if no points can be found in that area : stop the
       %computation
       
       if size(temp_neighours,1)==0
           keep_climb_gradient = false;
       end

        completed_iterations  = completed_iterations +  1 ;
    end
    List_modes{k,1} = mean_pos;
end
List_modes_table = [];


%Transforming the list into a table
for k=1:size(List_modes,1)
    List_modes_table = [List_modes_table ; List_modes{k}];
end

%%Now we need to aggregate the modes that are too close 
%%We can create some adjacency matrix and then convert it into a matrix
%The different connected components are then extracted

%First : range search
Adj_list = rangesearch(List_modes_table,List_modes_table,epsilon_distance);
Adj_list = cellfun(@(x) double(ismember(1:size(List_modes_table,1),x(:))),Adj_list,'UniformOutput',false);
Adj_list =cell2mat(Adj_list);

%Second : identification of the connected components
espilon_graph = graph(sparse(Adj_list));
Final_seeds_clustering = conncomp(espilon_graph);

%Third : geting the centroid of the meta clusters

Final_clusters = splitapply(@(x) median(x,1),List_modes_table,Final_seeds_clustering.');

%%Last steps : getting the clusters for all the points 

Final_clustering = Final_seeds_clustering(Attributed_seed);

end


