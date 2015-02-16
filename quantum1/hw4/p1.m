clear all; close all
N=100;
w0=1;
w = w0*2*sin( ((1:N).*(pi/2))./(N+1) );
k_w = [0.1:0.2:2]; %kT/hw0
E=zeros(1,length(k_w));
for kidx = 1:length(k_w)
    kwt = 1/(k_w(kidx));
    E(kidx)=sum( w./( exp(w.*kwt)-1 ) );
end

figure; plot(E/N);