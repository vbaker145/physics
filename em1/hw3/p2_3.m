clear all; close all

x = 0:0.01:5;
x0 = [2 1 1]; y0 = [1 1 2];
figure(10); hold on;
for jj=1:3   
    t1 = y0(jj)./((x-x0(jj)).^2+y0(jj)^2);
    t2 = y0(jj)./((x+x0(jj)).^2+y0(jj)^2);
    sigma_x = -(1/pi).*(t1-t2);
    subplot(1,3,jj); plot(x,sigma_x)
    xlabel('X','FontSize',12)
    ylabel('\sigma / \lambda','FontSize',12)
end