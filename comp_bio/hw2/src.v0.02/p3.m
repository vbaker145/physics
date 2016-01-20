clear all; close all

%Initial speed distribution
vi = load('out_initial_v.dat');
N = length(vi);
nbins = length(vi)/30;
[idist xi] = hist(vi, nbins);
idist = idist./N;
figure(20); subplot(1,2,1); bar(xi, idist); title('Initial speed distribution')
%Obtain Maxwell distribution 
%T=0.5 for 3D, T=0.7502
mdi = md(xi, 0.8, 0.7502 ,2);
hold on; plot(xi,mdi)
erri = mean((idist-mdi).^2)

%Final speed distribution
vf = load('out_final_v.dat');
[fdist xf] = hist(vf, nbins);
fdist = fdist./N;
figure(20); subplot(1,2,2); bar(xf, fdist); title('Final speed distribution')
%Obtain Maxwell distribution 
mdf = md(xf, 0.8, 0.5675 ,2);
hold on; plot(xf,mdf)
errf = mean((fdist-mdf).^2)