clear

n = 101;
x = linspace(-4,4,n);
delx2 = (x(2)-x(1))^2;
a = zeros(n,n);

%First row
a(1,1) = -2/delx2 - x(1)^2;
a(1,2) = 1/delx2;
for jj = 2:n-1
   a(jj, jj-1) = 1/delx2;
   a(jj, jj+1) = 1/delx2;
   a(jj, jj)   = -2/delx2-x(jj)^2;
end
%Last row
a(n,n-1) = 1/delx2;
a(n,n) = -2/delx2-x(n)^2;

[eVecs, eVals] = eig(a);

eVals = diag(eVals);

figure; plot(x, -eVecs(:,end-5:end));