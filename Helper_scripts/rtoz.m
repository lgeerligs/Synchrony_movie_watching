function [z]=rtoz(var)
z=zeros(size(var)).*NaN;
nonan=find(~isnan(var)&var~=0);
z(nonan)=0.5.*log((1+var(nonan))./(1-var(nonan)));
