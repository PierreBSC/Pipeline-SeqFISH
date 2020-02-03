function save_parameters_json(parameters,file_name)

% Create json string
jsonStr = jsonencode(parameters);

% Create folder if not present
[path, name] = fileparts(file_name);
if ~isfolder(path)
    mkdir(path)
end

% Save to file
fid = fopen(file_name, 'w');

if fid == -1
    error('Cannot create JSON file');
end
fwrite(fid, jsonStr, 'char');
fclose(fid);

disp('Parameters saved!')