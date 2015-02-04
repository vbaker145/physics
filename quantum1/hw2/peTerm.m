function [ d ] = peTerm( idx1, idx2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x = linspace(-4,4,1001);
sigma = abs(x(2)-x(1));
d = (sqrt(pi)/2)*sigma*exp(-(idx1-idx2)^2/4)*(((x(idx1)+x(idx2))/2)^2+0.5*sigma^2);

end

