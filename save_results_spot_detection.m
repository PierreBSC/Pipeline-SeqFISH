function save_results_spot_detection(Parameters,Analysis_result,path_save_base)

% Save results
for R=1:Parameters.N_round
    name_round = strcat("Round_",num2str(R));
    %mkdir(name_round_dir)
    for P=1:Parameters.N_position
        name_position = strcat("Position_",num2str(P));
        path_save = fullfile(path_save_base,name_round,name_position);
        mkdir(path_save)
        
        X = Analysis_result.Spot_analysis_raw{R,1,P};
        writematrix(X,fullfile(path_save,"Spot_position.txt"),"delimiter","\t");
    end
end
disp('Spot detection results saved')