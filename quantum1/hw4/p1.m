
N=100;
dt = 0.01; %Temperature spacing
k_w = (dt:dt:1); %kT/hw0
t=0:0.1:3; %Temperature
E=zeros(1,length(k_w));

w = 1*2*sin( ((1:N).*(pi/2))./(N+1) );    
for kidx = 1:length(k_w)
    E(kidx)=sum( w./(exp(w./k_w(kidx))-1) );
end
E=E/N; %Normalize to number of particles

%Calculate the specific heat with numerical derivative
sheat = diff(E)./0.01;
