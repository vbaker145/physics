clear all; close all
ys = [5 6 7 8 9 10 12 15 20];

for jj=1:10
   therm_name = ['src.v0.09/results/thermo_' num2str(jj) '.dat'];
   
   thermo = load(therm_name);
   r2 = thermo(end/2:end ,8);
   op(jj) = mean(r2);
   op_std(jj) = std(r2);
end

figure; errorbar(1:10, op, op_std)
