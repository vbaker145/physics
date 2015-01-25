%QM homework 1

alpha = [.9*.9]';
n = 1:100;

alpha_m = repmat(alpha,1,length(n));
n_m     = repmat(n, length(alpha), 1);
eff = (alpha_m.^n_m).*(cos(pi./(2*n_m)).^2).^n_m;
figure; plot(eff'); 
c = num2str(alpha); legend(c)