clear all; close all;

%TM modes
z = [2.41, 3.83, 5.14, 5.52];
p = [0 1 2];
rl = 0:0.01:2;
figure(10); subplot(1,2,1); hold on; title('TM modes'); xlabel('R/L')
for jj=1:length(p)
    for kk=1:length(z);
        omega(jj,kk,:) = sqrt(z(kk)^2+(p(jj)*pi.*rl).^2);
        plot(rl, squeeze(omega(jj,kk,:)));
    end
end

%TE modes
z = [1.84,3.05, 3.83, 5.33];
p = [1 2 3];
rl = 0:0.01:2;
subplot(1,2,2);  hold on; title('TE modes'); xlabel('R/L')
for jj=1:length(p)
    for kk=1:length(z);
        omega(jj,kk,:) = sqrt(z(kk)^2+(p(jj)*pi.*rl).^2);
        plot(rl, squeeze(omega(jj,kk,:)));
    end
end