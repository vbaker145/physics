clear all; close all

x = -10:0.1:10;
k = -10:0.1:10;

%Part A
u = exp(-abs(x)/2);
A = (1/2*pi).*(1./(1/4+k.^2) ).^2;
figure; subplot(2,1,1); plot(x,u); xlabel('x'); ylabel('|u(x,0)|^2');
subplot(2,1,2); plot(k,A); xlabel('k'); ylabel('|A(k)|^2');

%Part B
u = exp(-abs(x).^2/4);
A = 2*exp(-2.*k.^2);
figure; subplot(2,1,1); plot(x,u); xlabel('x'); ylabel('|u(x,0)|^2');
subplot(2,1,2); plot(k,A); xlabel('k'); ylabel('|A(k)|^2');

%Part C
xin = find(abs(x)<1);
xr = x(xin);
u = zeros(1,length(x));
u(xin) = 1-abs(xr);
A = (1/2*pi).*((1-cos(k))./k.^2).^2;
figure; subplot(2,1,1); plot(x,u); xlabel('x'); ylabel('|u(x,0)|^2');
subplot(2,1,2); plot(k,A); xlabel('k'); ylabel('|A(k)|^2');

%Part D
xin = find(abs(x)<1);
xr = x(xin);
u = zeros(1,length(x));
u(xin) = 1;
A = (2/pi).*(sin(k)./k).^2;
figure; subplot(2,1,1); plot(x,u); xlabel('x'); ylabel('|u(x,0)|^2');
subplot(2,1,2); plot(k,A); xlabel('k'); ylabel('|A(k)|^2');

