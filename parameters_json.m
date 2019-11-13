%% Create structure with parameters
parameters = Define_parameters();

%% Write to json
jsonStr = jsonencode(parameters);
fname= 'seqfish_parameters.json';
fid = fopen(fname, 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, jsonStr, 'char');
fclose(fid);

%% Read from json
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
parameters_read = jsondecode(str);

