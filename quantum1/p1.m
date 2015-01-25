clear

x = linspace(-4,4,1000);
h(:,1) = ones(1,length(x));
h(:,2) = 2*x;
h(:,3) = 4*x.^2-2;
h(:,4) = 8*x.^3-12*x;
h(:,5) = 16*x.^4-48*x.^2+12;
h(:,6) = 32*x.^5-160*x.^3+120*x;

for jj=1:6
   h(:,jj) = h(:,jj).*exp(-x.^2)';
   norm = 1/(2^(jj-1)*factorial(jj-1)*sqrt(pi));
   h(:,jj) = h(:,jj)*norm;
end

figure; plot(x, h);
