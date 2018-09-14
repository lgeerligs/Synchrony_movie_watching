function [r]=ztor(z)
% Inverse of the Fisher's transformation of the correlation coefficient
% The input z is the Fisher transformed correlation coefficient 
% The output r is the correlation coefficient
% L Geerligs 14-09-2018

r=zeros(size(z));

for i=1:numel(z)
    r(i)=(exp(2*z(i))-1)./(exp(2*z(i))+1);
end
   