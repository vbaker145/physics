clear all; close all
t=linspace(0,4*pi,100);
e1 = exp(-1i*(t));
e2 = exp(1i*(t));
p10 = (abs(e1+e2).^2)./4;
p01 = (abs(e1-e2).^2)./4;
figure(10); plot(t, p10);
hold on; plot(t, p01, 'r');
legend('P10', 'P01')
xlabel('time')
ylabel('probability')



