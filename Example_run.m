[Analysis_result,Parameters] = Create_experiment();

Analysis_result=Spot_detection(Analysis_result,Parameters);

Analysis_result=Compute_global_stitching(Analysis_result,Parameters,1,2);

[Analysis_result] = Spot_based_segmentation(Analysis_result,Parameters);

Spot_visualisation(Analysis_result,Parameters,1,1,4);

hold on
gscatter(Final_table(:,1),Final_table(:,2),Final_table(:,4))


X = cellfun(@(x) size(x,1) ,Analysis_result.Spot_analysis_raw)
X = squeeze(X(:,1,:));



[Analysis_result,Parameters] = Create_experiment();
Analysis_result=Spot_detection(Analysis_result,Parameters);
Spot_visualisation(Analysis_result,Parameters,10,2,1);

cd
for R=1:Parameters.N_round
    name_round_dir = strcat("Round_",num2str(R));
    mkdir(name_round_dir)
    for P=1:Parameters.N_position
        name_position = strcat(name_round_dir,"/Position_",num2str(P),"/");
        mkdir(name_position)

        X = Analysis_result.Spot_analysis_raw{R,1,P};
        writematrix(X,strcat(name_position,"Spot_position.txt"),"delimiter","\t");
    end
end
