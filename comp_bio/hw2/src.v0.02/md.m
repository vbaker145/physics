function [ md ] = md( vs, rho, T, d )
%Obtain Maxwell distribution 
md = 4*pi*rho*(1/(2*pi*T))^(d/2).*(vs.^(d-1)).*exp(-(vs.^2)./(2*T));
%Normalize?
md = md./sum(md);

end

