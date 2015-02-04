function [ d ] = keTerm(idx1, idx2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x = linspace(-4,4,1001);
sigma = abs(x(2)-x(1));
d = (sqrt(pi)/sigma)*exp(-(idx1-idx2)^2/4)*(-( (idx1+idx2)^2 )/4 ...
                              +0.5*sigma^2+x(idx1)*x(idx2));

end

