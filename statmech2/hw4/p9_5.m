clear all; close all;

V = [100];
z=0:0.01:5;
alpha = 3;
Vm = repmat(V',1, length(z));
zm = repmat(z,length(V), 1);
p1 = log(1+zm);
p2 = (1./Vm).*log(1+zm.^(Vm));
p_kt=log(1+zm)+(1./Vm).*log(1+zm.^(Vm.*alpha));

figure; plot(z, p_kt')
xlabel('z'); ylabel('P(z)')

lv = (zm./(zm+1))+(alpha./(zm.^(-alpha.*Vm)+1));
figure; plot(z,lv');
ylim([0 5]);
xlabel('z'); ylabel('v(z)')
