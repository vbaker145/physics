clear all; close all
%g845 = load('results/gr_845.dat');
g845 = load('results/gr_845_3D.dat');

gr = g845(:,1);
g = g845(:,2);
p = 0.844;

sidx = find(gr==min(gr));
sidx = sidx(end); %Only use last data set (guaranteed equilibrium)
g = g(sidx:end);
gd = gr(sidx:end);

%Assume highest peak is NN distance
[mv maxIdx] = max(g);
gp = g(maxIdx:end);
%Find next minima
gp(1:end-1) = gp(1:end-1)-gp(2:end); 
idx = find(gp<0);
minIdx = maxIdx+min(idx);
figure; plot(gd, g,'x-'); title(['G (p=' num2str(p) ')']);
hold on; plot(gd(maxIdx),g(maxIdx),'go','MarkerSize',10);
hold on; plot(gd(minIdx),g(minIdx),'ro','MarkerSize',10);
dr = diff(gd(1:minIdx+1));

%2-D integral
%gint(jj) = p(jj)*2*pi*(sum(g(1:minIdx).*diff(gd(1:minIdx+1)).*gd(1:minIdx)));

%3-D integral 
gint = p*4*pi*(sum(g(1:minIdx).*diff(gd(1:minIdx+1)).*gd(1:minIdx).^2));

%3-D integral up to L=2
drMean = mean(diff(gd));
gtot = p*4*pi*drMean*(sum(g.*gd.^2));

guni = (4/3)*pi*p*(7.5/2)^3

gbox = p*7.5^3



