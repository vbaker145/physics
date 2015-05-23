clear all; close all;

%Construct the matrix
N=1000; %Must be even
Nruns = 50;
e1 = []; er=[];
for runIdx = 1:Nruns
    M = zeros(N,N);
    for jj=2:N
        for kk=1:jj-1
            M(jj,kk) = rand()-1;
            M(kk,jj) = M(jj,kk);
        end
    end
    ev = sort(eig(M));
    e1(runIdx) = ev(1);
    er(runIdx,:) = ev(2:end);
    disp(['Run ' num2str(runIdx)])
end

figure; plot(e1,'rx')
hold on; plot(min(er'),'bx')
hold on; plot(mean(er'),'gx')
hold on; plot(max(er'),'mx')
xlabel('Trial #')
ylabel('Eigenvalue')
legend('Minimum eigenvalues', 'Next-lowest eigenvalues', ...
            'Mean eigenvalues', 'Max eigenvalues');
