function [Concat_Image] = LoadImage(input_dir,load_specific_file,specific_file)
%LoadImage : loading and aggregation of the images 
%   Each image coming from a specific round and zone is encoded by 4/5
%   different tiff images, one by chanel.

if nargin<3
    load_specific_file = false;
    specific_file = [];
end

list_files=dir([input_dir, '/*.tif']); %%Only taking tif files to  read

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
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;


Concat_Image=zeros(nImage,mImage,numstack,l,'single');

for i=1:l
    for j=1:numstack 
        temp_im=im2single(imread(fnames{i},j)); %%% We read image by image
        Concat_Image(:,:,j,i)=temp_im;
    end
end

disp("Tiff files read and concatenated");

end

