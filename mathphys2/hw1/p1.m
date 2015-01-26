clear all

n = 32768; 
theta = 2*pi*(0:n-1)*(1/n); 
r = rand(n,1)*2-1;
r = r.*hamming(n);
dss = 0;
theta_s = sum(cos(theta).^2);
for jj=1:n
    dss = dss + cos(theta(jj))^2.*theta_s;
end
dsc = sum(cos(theta).^4);

