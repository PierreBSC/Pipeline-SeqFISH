function save_parameters_json(parameters,file_name)

% Save detection settings
jsonStr = jsonencode(parameters);

fid = fopen(file_name, 'w');

if fid == -1
    error('Cannot create JSON file');
end
fwrite(fid, jsonStr, 'char');
fclose(fid);

disp('Parameters saved!')