function [r]=ztor(z)
r=zeros(size(z));

for i=1:numel(z)
    r(i)=(exp(2*z(i))-1)./(exp(2*z(i))+1);
end
   