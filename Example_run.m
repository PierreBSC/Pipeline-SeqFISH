[Analysis_result,Parameters] = Create_experiment()
Analysis_result=Spot_detection(Analysis_result,Parameters)
[Analysis_result] = Spot_based_segmentation(Analysis_result,Parameters,true)

Spot_visualisation(Analysis_result,Parameters,1,1,1)
Segmentation_visualisation(Analysis_result,Parameters,1,1,1,true)