clear all; close all

%Define constants
hbar = 6.626e-34/(2*pi);
k = 1.3806e-23;
t = 1e-3:1:300; %Plot from 0 to 300 Kelvin
beta = 1./(k.*t);

omega = 2*pi*4.2e9; %Microwave frequency
n_m = 1./(exp(beta*hbar*omega)-1);

omega = 2*pi*500e12;
n_o = 1./(exp(beta*hbar*omega)-1);

figure(1);
subplot(2,1,1); plot(t,n_m); 
xlabel('Temperature (K)','FontSize', 12);
ylabel('n_{e\omega}','FontSize', 14)
subplot(2,1,2); plot(t,n_o);
xlabel('Temperature (K)','FontSize', 12)
ylabel('n_{eo}','FontSize', 14)

