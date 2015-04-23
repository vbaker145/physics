clear alll; close all

%Input 1970 particle data
eDelta = [1231 1233 1235 1237];
eDeltaStd = [10 10 10 10];
eSigma = [1381 1385 1389];
eSigmaStd = [11 11 11];
eRho = [1530 1540];
eRhoStd = [15 15];

d1 = mean(eSigma)-mean(eDelta);
d2 = mean(eRho)-mean(eSigma);

eOmega = mean(eRho)+mean(d1,d2);

%3D harmonic oscillator potentials, N=3
states= [3 0 0; 2 1 0; 2 0 1; 1 2 0; 1 0 2; 1 1 1; 0 3 0; 1 2 0; 0 2 1; 0 0 3];

bm = 257.33;
m = [(bm+150) bm bm];
e = (states(:,1)+1/2)*m(1) + (states(:,2)+1/2)*m(2)+(states(:,3)+1/2)*m(3);
%figure; plot(e,'x');

npert = 50000; maxPert = 200;
pert = (rand(npert, 3)-1/2) * maxPert;
mp = repmat(m,npert,1) + pert;
%pert = -10:0.1:10;
%mp = repmat(m,length(pert)*3 ,1);
%mp(1:length(pert), 1) = m(1)+pert;
%mp(length(pert)+1:length(pert)*2,2) = m(2)+pert;
%mp(length(pert)*2+1:length(pert)*3,3) = m(3)+pert;

ep = []; mmse = [];
[actVals, actIdx] = sort([eDelta eSigma eRho]);
actVar = [eDeltaStd eSigmaStd eRhoStd ].^2;
actVar = actVar(actIdx);

for jj=1:length(mp)
  ep(jj,:) = sort( (states(:,1)+1/2)*mp(jj,1) + (states(:,2)+1/2)*mp(jj,2)+(states(:,3)+1/2)*mp(jj,3) );
  mmse(jj) = mean((ep(jj,1:9)-actVals).^2);
end

[mv idx] = min(mmse);
mv
figure; plot(ep(idx,:),'x')
hold on; plot(actVals, 'rx')

best = mp(idx,:);
best10 = (states(:,1)+1/2)*best(1) + (states(:,2)+1/2)*best(2)+(states(:,3)+1/2)*best(3);
chi2 = 0;
for jj=1:9
    chi2 = chi2 + (ep(idx,jj)-actVals(jj))^2/actVar(jj);
end

chi2_no4 = 0;
for jj=[1 2 3 5 6 7 8 9]
    chi2_no4 = chi2_no4 + (ep(idx,jj)-actVals(jj))^2/actVar(jj);
end

