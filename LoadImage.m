function [Concat_Image] = LoadImage(input_dir,load_specific_file,specific_file)
%LoadImage : loading and aggregation of the images 
%   Each image coming from a specific round and zone is encoded by 4/5
%   different tiff images, one by chanel.

if nargin<3
    load_specific_file = false;
    specific_file = [];
end

list_files=dir([input_dir, '/*.tif']); %%Only taking tif files to  read

if size(list_files,1)==0
    list_files=dir([input_dir, '/*.tiff']); %%Only taking tif files to  read
end


if load_specific_file
    list_files = list_files(specific_file);
end

l = size(list_files);
l = l(1);
fnames = cell(l,1); %% 

for j=1:l
        fnames{j} = [input_dir,list_files(j).name]; 
end

numstack = size(imfinfo(fnames{1}),1);

InfoImage=imfinfo(fnames{1});
Bit_depth = InfoImage(1).BitDepth;
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;


Concat_Image=zeros(nImage,mImage,numstack,l,'double');

for i=1:l
    for j=1:numstack 
        temp_im=(imread(fnames{i},j)); %%% We read stack by stack
        temp_im = im2double(int16(temp_im));

        Concat_Image(:,:,j,i) = (temp_im);
    end
end

disp("Tiff files read and concatenated");

end

