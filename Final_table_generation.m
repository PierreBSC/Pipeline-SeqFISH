function [Analysis_result] =Final_table_generation(Analysis_result,Parameters)
%Final function for image processing : generate the two tables : one with
%cell position and one with cell gene expression.

for P=1:Parameters.N_position
    
    Cell_segmentation = Analysis_result.Spot_based_segmentation{P};
    N_cells = size(unique(Cell_segmentation(:,3)),1);
    
    for k=1:N_cells
        Cell_localisation = Cell_segmentation(Cell_segmentation(:,3)==k,1:2);
    end
    
end

end

