function output = double_dot(tensor1,tensor2)
%double_dot return the double_dot product of two tensors

if size(tensor1,1)~=size(tensor2,1)
    warning('dimensions are different');
else 
    output=sum(dot(tensor1,tensor2));   
end
end

