clear all;

N=9;

%Useful terms for calculating Fourier coefficients
Ag = [2*pi 0.2473 0.1557 0.0483 -0.0322];
Ag = [Ag(end:-1:2) Ag];
Bg = [0.2473 3.2194 0.1478 0.0617 -0.0066];
Bg = [Bg(end:-1:2) Bg];

%Create lattice vector (2D->1D)
[x, y] = meshgrid(-(N-1)/2:(N-1)/2, -(N-1)/2:(N-1)/2);
xidx = x+(N+1)/2; yidx = y+(N+1)/2; %Create matrix indices
Vg = 3*Ag(yidx)*Ag(xidx)+2*Ag(yidx)*Bg(xidx)+2*Ag(xidx)*Bg(yidx)+Bg(xidx)*Bg(yidx);
Vg = Vg./(4*pi^2);
%Create G matrices
Gx = 2*pi*x; Gy = 2*pi*y;                   
%Flatten matrices to vectors (MATLAB ordering)
Vg = reshape(Vg,1,[]);
Gx = reshape(Gx,1,[]);
Gy = reshape(Gy,1,[]);
%Create coefficient matrix
Vmat=[];
for jj=1:length(Vg);
   Vmat(jj,:) = circshift(Vg, jj-1, 2 ); 
end

%First section, (0,0)->(0,pi)
ky=0:pi/100:pi;
kx=zeros(1,length(ky));
ev1=[];
for jj=1:length(ky)
    Kmat = ((kx(jj)+Gx).^2+(ky(jj)+Gy).^2);
    Kmat = diag(Kmat)+Vmat;
    ev1(jj,:) = sort(abs(eig(Kmat)));
end
%figure; plot(ky, abs(ev1(:,1:4)))

%Second section (0,pi)->(pi,pi)
kx=0:pi/100:pi;
ky=ones(1,length(kx))*pi;
ev2=[];
for jj=1:length(ky)
    Kmat = ((kx(jj)+Gx).^2+(ky(jj)+Gy).^2);
    Kmat = diag(Kmat)+Vmat;
    ev2(jj,:) = sort(abs(eig(Kmat)));
end
%figure; plot(kx, abs(ev2(:,1:4)))

%Third section (0,pi)->(pi,pi)
kx=pi:-pi/100:0;
ky=kx;
ev3=[];
for jj=1:length(ky)
    Kmat = ((kx(jj)+Gx).^2+(ky(jj)+Gy).^2);
    Kmat = diag(Kmat)+Vmat;
    ev3(jj,:) = sort(abs(eig(Kmat)));
end
%figure; plot(abs(ev3(:,1:4)))

evt = [ev1(:,1:4); ev2(:,1:4); ev3(:,1:4)];
figure; plot(evt)
