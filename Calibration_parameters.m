%%Function to test the effect of changes on parameter values

[Analysis_result,Parameters] = Create_experiment();

Parameters.sigma_max = 5
Parameters.sigma_small = 1
S_range = 1:15;
h_meta_range = linspace(0.1,1.5,15);

Matrix_n_points = zeros(15);
Matrix_n_points_filtered = zeros(15);


writetable(array2table(Matrix_n_points),'/home/pbost/Desktop/Calibration_seqFISH/Matrix_n_points.txt','delimiter','\t')
writetable(array2table(Matrix_n_points_filtered),'/home/pbost/Desktop/Calibration_seqFISH/Matrix_n_points_filtered.txt','delimiter','\t')

Values_calibration = [S_range.',h_meta_range.'];
writetable(array2table(Values_calibration),'/home/pbost/Desktop/Calibration_seqFISH/Values_calibration.txt','delimiter','\t')




%%%Conclusion : global balance between S value and h_meta parameter... 
%%% S ~= 15*h_meta
%%%The use of spatial filtering : reduce the variations 


Parameters.S = 10;
Parameters.h_meta_parameter = 0.35;
N_pixels_sampled_range = round(linspace(5000,100000,15));

N_points_total =zeros(15,1);
N_points_total_filtered =zeros(15,1);

for i=1:size(N_pixels_sampled_range,2)
        Parameters.N_points_sample = N_pixels_sampled_range(i);
        disp(strcat(num2str(i)," ",num2str(j)))
        Analysis_result=Spot_detection(Analysis_result,Parameters);
        N_points_total(i) = size(Analysis_result.Spot_analysis_raw{4},1);
        N_points_total_filtered(i) = size(Analysis_result.Spot_analysis_filtered{4},1);
end

writematrix([N_pixels_sampled_range.',N_points_total,N_points_total_filtered],'/home/pbost/Desktop/Calibration_seqFISH/Sampling_saturation.txt','delimiter','\t')

plot(N_pixels_sampled_range,N_points_total)
hold on
plot(N_pixels_sampled_range,N_points_total_filtered)



