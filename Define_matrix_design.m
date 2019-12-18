function Matrix_design = Define_matrix_design()

prompt = {  'Number of rounds',
            'Number of channels'
};
dlgtitle = 'Experiment design information';
dims = [1 35];
definput = {'1' '4'};

answer = inputdlg(prompt,dlgtitle,dims,definput);




if ~ isempty(answer)
    
    % Read parameters
    N_Rounds = str2num(answer{1}); 
    N_Channels = str2num(answer{2});
    
    Table_base = repmat("RNA",1,N_Rounds*N_Channels); 
    Table_base=cellstr(Table_base);
    
    AddOpts.Resize='on';
    AddOpts.WindowStyle='normal';
    AddOpts.Interpreter='tex';

    prompt = 1:N_Rounds*N_Channels;
    prompt = num2cell(prompt);
    prompt = cellfun(@num2str,prompt,'UniformOutput',false);
    Matrix_design=inputdlgcol(prompt,'Provide Matrix Design',[1,34],Table_base,AddOpts,N_Channels);
    Matrix_design = reshape(Matrix_design,N_Rounds,N_Channels);
    Matrix_design =  cell2table(Matrix_design);
    
    names_col = cell(0);
    for i=1:N_Channels
        names_col = [names_col,strcat('Channel_',num2str(i))];
    end
    Matrix_design.Properties.VariableNames = names_col;
    
else
    disp('Parameter definition canceled')
    Parameters = {};
end
    

    
end
