function [Final_Image] = Optimised_projection(Image,score_name,size_window)
%For a three-dimensional array : compute locally the best stack using a
%given measure of focus

if nargin < 3
    score_name = 'GLVN';
    size_window = 16;
end

Split_image = mat2cell(Image,repmat(size_window,1,size(Image,1)/size_window),repmat(size_window,1,size(Image,2)/size_window),size(Image,3));

Local_score = cellfun(@(x) get_stack_score(x,score_name),Split_image,'UniformOutput',false);
Best_stack = cellfun(@(x) find(x==min(x)),Local_score,'UniformOutput',false);
Best_stack = cellfun(@(x) x(1),Best_stack,'UniformOutput',true);

Final_Image = cell(size(Split_image));

for i=1:size(Split_image,1)
    for j=1:size(Split_image,2)
        x = Split_image{i,j};
        Final_Image{i,j} = x(:,:,Best_stack(i,j));
    end
end

Final_Image = cell2mat(Final_Image);

end
%************************************************************************

function focus_score = get_stack_score(Image,score_name)
focus_score = [];

for k = 1:size(Image,3)
  focus_score = [focus_score fmeasure(Image(:,:,k),score_name)];
end

end

