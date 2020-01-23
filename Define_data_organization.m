function Parameters = Define_data_organization(Parameters)

prompt = {  'Name of experiment', ...
            'Number of positions', ...
            'Number of rounds', ...
            'Number of channels'
};
dlgtitle = 'Experiment design information';
dims = [1 35];

if isfield(Parameters,'Experiment_name')
    definput = {Parameters.Experiment_name, ...
                num2str(Parameters.N_position), ...
                num2str(Parameters.N_round), ...
                num2str(Parameters.N_channel), ...
                };
else
    definput = {'seq-FISH analysis','1','1','4'};
end

answer = inputdlg(prompt,dlgtitle,dims,definput);

if ~ isempty(answer)
    
    % Read parameters
    Experiment_name = answer{1}; 
    N_position = str2num(answer{2});
    N_Rounds = str2num(answer{3}); 
    N_Channels = str2num(answer{4});
    
    Table_base = repmat("RNA",1,N_Rounds*N_Channels); 
    Table_base=cellstr(Table_base);
    
    AddOpts.Resize='on';
    AddOpts.WindowStyle='normal';
    AddOpts.Interpreter='tex';

    prompt = 1:N_Rounds*N_Channels;
    prompt = num2cell(prompt);
    prompt = cellfun(@num2str,prompt,'UniformOutput',false);
    Matrix_design = inputdlgcol(prompt,'Provide Matrix Design',[1,34],Table_base,AddOpts,N_Channels);
    Matrix_design = reshape(Matrix_design,N_Rounds,N_Channels);
    Matrix_design =  cell2table(Matrix_design);
    
    names_col = cell(0);
    for i=1:N_Channels
        names_col = [names_col,strcat('Channel_',num2str(i))];
    end
    Matrix_design.Properties.VariableNames = names_col;
    
    Parameters.Matrix_design = Matrix_design;
    Parameters.Experiment_name = Experiment_name;
    Parameters.N_position = N_position;
    Parameters.N_round = size(Matrix_design,1);
    Parameters.N_channel = size(Matrix_design,2);
else
    disp('Definition of experiment canceled.')
end







