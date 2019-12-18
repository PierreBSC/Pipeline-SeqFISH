function Segmentation_design = Select_segmentation_genes(Parameters)


    
    % Read parameters
    N_Rounds = size(Parameters.Matrix_design,1);
    N_Channels = size(Parameters.Matrix_design,2);
    
    Table_base = repmat("Use",1,N_Rounds*N_Channels); 
    Table_base=cellstr(Table_base);
    
    AddOpts.Resize='on';
    AddOpts.WindowStyle='normal';
    AddOpts.Interpreter='tex';

    prompt = 1:N_Rounds*N_Channels;
    prompt = num2cell(prompt);
    prompt = cellfun(@num2str,prompt,'UniformOutput',false);
    Segmentation_design = inputdlgcol(prompt,'Provide which channels to use for clustering',[1,34],Table_base,AddOpts,N_Channels);
    Segmentation_design = reshape(Segmentation_design,N_Rounds,N_Channels);
    Segmentation_design =  cell2table(Segmentation_design);
    

end
