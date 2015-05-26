clear all;

N=9;
G=-(N-1)/2:(N-1)/2;
%G=[0:(N-1)/2 -(N-1)/2:-1]; %Alternate order of eigenvectors
G = G*2*pi;
k=0:pi/N:pi;

Vg = [2.8099  0.2685 0.1579 0.0486 -0.0323];
Vg = [Vg(5:-1:2) Vg];
%Vg(6:9) = Vg(5:-1:2); %Alternate order of eigenvectors
Vg = Vg*(1/(2*pi))*2.5;
for jj=1:N
   Vmat(jj,:) = circshift(Vg, jj-1, 2 ); 
end
ev = [];
for jj=1:length(k)
    Kmat = 0.5*(k(jj)+G).^2;
    Kmat = diag(Kmat)+Vmat;
    ev(jj,:) = sort(eig(Kmat));   
end

figure; plot(k, real(ev(:,1:4)))
