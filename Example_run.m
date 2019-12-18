[Analysis_result,Parameters] = Create_experiment();
Analysis_result=Spot_detection(Analysis_result,Parameters);
[Analysis_result] = Spot_based_segmentation(Analysis_result,Parameters);

Spot_visualisation(Analysis_result,Parameters,1,1,2)
Segmentation_visualisation(Analysis_result,Parameters,1,1,4,true)

Compute_spots_QC_plot(Analysis_result,Parameters)

[Analysis_result,Concatenated_RNA_data] = Compute_cell_RNA_expression(Analysis_result,Parameters);
[Analysis_result,Concatenated_Flu_data] = Compute_cell_fluorescence(Analysis_result,Parameters,'IF');
[General_cell_info] = Compute_cell_properties(Analysis_result,Parameters);


Table_Cy5 = [];

for P=1:Parameters.N_position
    X = Analysis_result.Spot_analysis_raw{1,1,P};
    X(:,3) = X(:,3)*5;
    P_rep = repmat(P,size(X,1),1);
    X = [X , P_rep];
    Table_Cy5 = [Table_Cy5;X];
end

Table_Cy5 = table(Table_Cy5);

writetable(Table_Cy5,'/home/pbost/Desktop/Hela_cropped_Cy5.txt','delimiter','\t')



Table_Cy3 = [];

for P=1:Parameters.N_position
    X = Analysis_result.Spot_analysis_raw{2,1,P};
    X(:,3) = X(:,3)*5;
    P_rep = repmat(P,size(X,1),1);
    X = [X , P_rep];
    Table_Cy3 = [Table_Cy3;X];
end

Table_Cy3 = table(Table_Cy3);

writetable(Table_Cy,'/home/pbost/Desktop/Hela_cropped_Cy3.txt','delimiter','\t')