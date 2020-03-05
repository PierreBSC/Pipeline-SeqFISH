function [Filtered_spots] = Spot_filtering(Unfiltered_spots,Spot_intensity,Intensity_threshold,Autolag_score)
%Filtering of the spot based on various criteria

l = size(Unfiltered_spots,1);
Filtered_spots = cell(1,l);
for k=1:l
    
    temp_intensity = Spot_intensity{k};
    
    selected_spots = temp_intensity.Variables > Intensity_threshold(k);
    
    temp_position = Unfiltered_spots{k};
    temp_position = temp_position(selected_spots,:);
        
	Filtered_spots{k} = temp_position;

end

end
