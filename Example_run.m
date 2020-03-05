[Analysis_result,Parameters] = Create_experiment();

Parameters.N_scales = 3;
Parameters.sigma_small = 0.8;
Parameters.sigma_max = 2;
Parameters.Spot_detection_method = "Multiscale";
Parameters.Quantile_parameter = 0.00001;
Parameters.background_sigma_parameter = 100;

Analysis_result=Spot_detection(Analysis_result,Parameters);
Analysis_result=Compute_global_stitching(Analysis_result,Parameters,1,2);
[Global_analysis_results] = Compute_global_point_position(Analysis_result,Parameters)

Parameters.N_position = 1;

[Global_analysis_results] = Spot_based_segmentation(Global_analysis_results,Parameters);

Spot_visualisation(Analysis_result,Parameters,1,2,2);

hold on
gscatter(Final_table(:,1),Final_table(:,2),Final_table(:,4))


X = cellfun(@(x) size(x,1) ,Analysis_result.Spot_analysis_raw)
X = squeeze(X(:,1,:));


[Analysis_result,Parameters] = Create_experiment();
Analysis_result=Spot_detection(Analysis_result,Parameters);
Spot_visualisation(Analysis_result,Parameters,1,1,1);


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
