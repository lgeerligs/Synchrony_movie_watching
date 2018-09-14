function [z]=rtoz(r)
% Fisher's transformation of the correlation coefficient
% The input r is the correlation coefficient
% The output z is the Fisher transformed correlation coefficient 
% L Geerligs 14-09-2018

z=zeros(size(r)).*NaN;
nonan=find(~isnan(r)&r~=0);
z(nonan)=0.5.*log((1+r(nonan))./(1-r(nonan)));
