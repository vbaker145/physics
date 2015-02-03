clear; close all

x = linspace(-4,4,1001);
npts = length(x);
sigma = abs(x(2)-x(1));

%Form KE+PE matrix
KM = zeros(npts,npts);
PM  = zeros(npts,npts);
PHM = zeros(npts,npts);
EM  = zeros(npts, npts);
for ii=1:npts
   for jj=max(1,ii-5):min(npts,ii+5)
       %KE contribution
       KM(ii,jj) = (sqrt(pi)/sigma)*exp(-(ii-jj)^2/4)* ... 
                              (-(x(ii)+x(jj))^2/4 ...
                              +0.5*sigma^2+x(ii)*x(jj));
       %PE contribution
       PM(ii,jj) = (sqrt(pi)/2)*sigma*exp(-(ii-jj)^2/4)*(((x(ii)+x(jj))/2)^2 ...
                        +0.5*sigma^2);

       %Cross function             
       EM(ii,jj) = sqrt(pi)*sigma*exp(-(ii-jj)^2/4);
   end
end
KPM = KM+PM;
[eVecs, eVals] = eig(KPM, EM); %Generalized eigenvals
eVals = diag(eVals);
figure; plot(x, eVecs(:,1:6).^2)
