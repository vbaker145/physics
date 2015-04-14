clear all; close all
t=linspace(0,4*pi,100);
theta = atan(2);
e1 = (cos(theta/2)^2)*exp(-1i*(t));
e2 = (sin(theta/2)^2)*exp(1i*(t));

p10 = (e1+e2).*conj(e1+e2);
p01 = 1-p10;

figure(10); plot(t, p10);
hold on; plot(t, p01, 'r');
legend('P10', 'P01')
xlabel('time')
ylabel('probability')




