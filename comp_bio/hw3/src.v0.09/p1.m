clear all; close all
g845 = load('results/gr_845.dat');
%g845 = load('results/gr_845_3D.dat');
g8 = load('results/gr_8.dat');
g7 = load('results/gr_7.dat');
g6 = load('results/gr_6.dat');
g5 = load('results/gr_5.dat');
gr = [g845(:,1) g8(:,1) g7(:,1) g6(:,1) g5(:,1)]; 
gv = [g845(:,2) g8(:,2) g7(:,2) g6(:,2) g5(:,2)]; 

p = [0.845 0.8 0.7 0.6 0.5];

for jj = 1:length(p)
    g = gr(:,jj);
    sidx = find(g==min(g));
    sidx = sidx(end); %Only use last data set (guaranteed equilibrium)
    g = gv(sidx:end,jj);
    gd = gr(sidx:end,jj);
    
    %Assume highest peak is NN distance
    [mv maxIdx] = max(g);
    gp = g(maxIdx:end);
    %Find next minima
    gp(1:end-1) = gp(1:end-1)-gp(2:end); 
    idx = find(gp<0);
    minIdx = maxIdx+min(idx);
    figure; plot(gd, g,'x-'); title(['G (p=' num2str(p(jj)) ')']);
    hold on; plot(gd(maxIdx),g(maxIdx),'go','MarkerSize',10);
    hold on; plot(gd(minIdx),g(minIdx),'ro','MarkerSize',10);
    dr = diff(gd(1:minIdx+1));
    
    %2-D integral
    gint(jj) = p(jj)*2*pi*(sum(g(1:minIdx).*diff(gd(1:minIdx+1)).*gd(1:minIdx)));
    
    %3-D integral
    %gint(jj) = p(jj)*4*pi*(sum(g(1:minIdx).*diff(gd(1:minIdx+1)).*gd(1:minIdx).^2));
end

figure; plot(p, gint, 'x-')
xlabel('Density', 'FontSize', 14); ylabel('Coordination function','FontSize', 14);
set(gca,'FontSize',12)


