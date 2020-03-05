 
function [Final_clusters] = Diffusion_map(Spot_position,N_neighbors,N_comp,T)
%Diffusion map based cell segmentation
%   Function to segment  cells based on spot location : the approach is to the one used for
%   diffusion map analysis

if nargin < 4
    N_neighbors = 10; %%Number of neighbors used for the Mutual nearest neighbors graph
    N_comp = 50; %Number of eigenvalues to compute, reduce if too slow and increase if too many cells
    T = 15;
end

N_neighbors = min(N_neighbors,size(Spot_position,1));
Spot_position = Spot_position(:,1:2);
%%Local sigma is computed according to the paper "Self Tuning Spectral Clustering" 
%%See section "Local scaling"
[ ~ ,  knn_dist ]= knnsearch(Spot_position,Spot_position,'K',N_neighbors+1);
local_sigma = knn_dist(:,N_neighbors) ;
local_sigma = local_sigma * local_sigma.' ; 

Dist_matrix = pdist2(Spot_position,Spot_position);
%%L : affinity matrix
L = exp( - Dist_matrix.^2 ./(2*local_sigma));

%%D : degree/normalisation matrix

D = diag(sum(L,2));

%%Normalisation of the affinity matrix
M = inv(D) *  L ;

%%Computation of the eigen values
%Checking that we don't have more eigenvalues to compute than total spot
%number
N_comp = min(N_comp,size(Spot_position,1));


[psi lambda_vector]= eigs(sparse(double(M)),N_comp,'largestabs');

lambda_vector = diag(lambda_vector);
Lambda_change = lambda_vector(1:(size(lambda_vector,1)-1)) ./ lambda_vector(2:(size(lambda_vector,1)));

[~ , N_eigen_values ] = max(Lambda_change);
N_eigen_values = N_eigen_values ;
Psi_normalised = psi(:,1:N_eigen_values);
Psi_normalised = arrayfun(@(x) x/sqrt(sum(x.^2)),Psi_normalised);
Psi_normalised = Psi_normalised .* transpose(lambda_vector(1:N_eigen_values).^T);
Final_clusters = kmeans(Psi_normalised,N_eigen_values,'Replicates',50);
 end
