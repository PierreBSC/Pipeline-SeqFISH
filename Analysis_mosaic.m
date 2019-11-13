Images_1 = LoadImage('/home/pbost/Desktop/Zstacks/mosaic_T_cells_1/');
Images_1 = Pre_processing(Images_1);

Images_2 = LoadImage('/home/pbost/Desktop/Zstacks/mosaic_T_cells_2/');
Images_2 = Pre_processing(Images_2);

Images_3 = LoadImage('/home/pbost/Desktop/Zstacks/STAT1_only_1/');
Images_3 = Pre_processing(Images_3);

Images_4 = LoadImage('/home/pbost/Desktop/Zstacks/STAT1_only_2/');
Images_4 = Pre_processing(Images_4);

List_score = cell(4,1);
List_images = {Images_1,Images_2,Images_3,Images_4};


for k=1:4
    
    temp_Image = List_images{k};
    h = fspecial('log', [30 30], 10);
    
    STAT1 = median(temp_Image(:,:,:,1),3);
    mCherry = median(temp_Image(:,:,:,2),3);
    
    STAT1_filter = imfilter((STAT1),h);
    STAT1_filter = -STAT1_filter;
    STAT1_filter = imadjust(STAT1_filter);
    
    h_parameter = 0.4;

    STAT1_Hdome =  imreconstruct(STAT1_filter-h_parameter,STAT1_filter,4);
    STAT1_Hdome = STAT1_filter - STAT1_Hdome;
        
    STAT1_enveloppe = Nuclei_identification(STAT1_Hdome);

    figure, imshow(imadjust(STAT1))
    hold on

    for i=1:size(STAT1_enveloppe) 
        x =polyshape(STAT1_enveloppe{i});
        plot(x);
    end
    hold off
    
    temp_score_list = [];
    
    for i=1:size(STAT1_enveloppe)
    
        temp_enveloppe = STAT1_enveloppe{i};
        temp_enveloppe_size = area(polyshape(temp_enveloppe));

        temp_enveloppe_ROI = poly2mask(temp_enveloppe(:,1),temp_enveloppe(:,2),size(STAT1,1),size(STAT1,2));
        global_STAT_signal = STAT1(temp_enveloppe_ROI);
        temp_nuclei_signal = mCherry(temp_enveloppe_ROI);
        
        temp_score = corr(global_STAT_signal,temp_nuclei_signal,"Type","Pearson");

        temp_score_list = [temp_score_list temp_score];

    end
    List_score{k} = temp_score_list;


end

merged_overlapping_score = [List_score{1} List_score{2} List_score{3} List_score{4}];

score_condition = [repmat("Mosaic_Tcell_1",size(List_score{1},2),1);repmat("Mosaic_Tcell_2",size(List_score{2},2),1);repmat("Control_1",size(List_score{3},2),1);repmat("Control_2",size(List_score{4},2),1)];

close all
boxplot(merged_overlapping_score,score_condition);

