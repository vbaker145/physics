clear all; close all

l = 1:.01:10;

f1 = -1./(l.^3);

f2 = 2*(1./(sqrt(l.^2+1))-1./l);

figure; plot(l,f1)
hold on; plot(l,f2,'r')
xlabel('x-y distance (a)','FontSize',12)
ylabel('Potential')
set(gca,'FontSize',12)