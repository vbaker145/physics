%Sketch electric field energy density
clear all; close all
dl = 0.0001;
l = 1:dl:2;

figure(10); hold on;
%Plates
Ep = ones(1,length(l));
subplot(1,3,1); plot(l,Ep); title('Plates')

%Spheres
Es = 1./(l.^4);
subplot(1,3,2); plot(l,Es); title('Spheres')

%Cylinders
Ec = 1./(l.^2);
subplot(1,3,3); plot(l,Ec); title('Cylinders')