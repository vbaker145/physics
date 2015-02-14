clear all; close all
N=100;

w = 2*sin( ((1:N).*(pi/2))./(N+1) );
k_w = 1./[0.1 0.5 1 10];
t=(0:0.01:1); %Actuall 1/t since exp term is 1/kT
E=zeros(length(t), length(k_w));
for kidx = 1:length(k_w)
    for tidx=1:length(t)
        E(tidx, kidx)=sum( (w.*k_w(kidx))./(exp(w.* (k_w(kidx)*t(tidx))) -1) );
    end
end

figure; plot(E/N);