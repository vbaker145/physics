clear all; close all

N=1000000;
u1 = rand(N,1);
[hu xu] = hist(u1,100);
g1 = zeros(N,1);
g1(1:2:end) = sqrt(-2.*log(u1(1:2:end))).*cos(2*pi.*u1(2:2:end));
g1(2:2:end) = sqrt(-2.*log(u1(1:2:end))).*sin(2*pi.*u1(2:2:end));
[hg xg] = hist(g1,100);
figure(10); subplot(1,3,1);  bar(xu, hu./N); title('Uniform distribution')
figure(10); subplot(1,3,2); bar(xg, hg./N); title('Gaussian distribution (method 1)')

sigma = 0.5; 
avg = 3;
ga = g1.*sigma + avg;
[ha xa] = hist(ga,100);
figure(10); subplot(1,3,3); bar(xa, ha./N); title('Gaussian distribution (\sigma 0.5, mean 3.0)')

nuni = [6 8 10 12 14];
hp = -3:0.1:3;
for jj=1:length(nuni)
   unif(:,jj) = rand(N,1);
   rnc = conv(unif(:,jj),ones(nuni(jj),1));
   rnc = rnc(nuni(jj):end);
   normed(:,jj) = rnc-nuni(jj)/2;
   dist(:,jj) = hist(normed(:,jj),hp);
   dist(:,jj) = dist(:,jj)./(length(rnc)*0.1);
end

gauss = 1/(sqrt(2*pi))*exp((-hp.^2)./2);
figure; plot(hp, dist,'x-')
hold on; plot(hp, gauss, 'r');
legend('6', '8', '10', '12', '14', 'Gaussian');
    