clear all; close all
ys = [5 6 7 8 9 10 12 15 20];

for jj=1:length(ys)
   r2t_name = ['src.v0.09/results/r2t_' num2str(ys(jj)) '.dat'];
   therm_name = ['src.v0.09/results/thermo_' num2str(ys(jj)) '.dat'];
   
   r2t = load(r2t_name);
   thermo = load(therm_name);
   c = polyfit(r2t(:,1),r2t(:,2),1);
   d(jj) = c(1);
   avg_temp(jj) = mean(thermo(:,3));
end

figure; plotyy(ys, d, ys, avg_temp)

figure; plot(ys, d./6,'x-'); hold on; plot(ys, 1./ys,'rx-')