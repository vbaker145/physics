clear

x = linspace(-4,4,1001);
npts = length(x);
sigma = abs(x(2)-x(1));

exp_vec = exp(-((-5:5).^2)/4);

%Form KE+PE matrix
KPM = zeros(npts,npts);
EM  = zeros(npts, npts
for ii=1:npts
   for jj=max(1,ii-5):min(npts,ii+5)
       ev = exp_vec(ii-jj+6);
       %KE contribution
       KPM(ii,jj) = KPM(ii,jj)+(sqrt(pi)/sigma)*ev*(-((x(ii)+x(jj))/2)^2 ...
                              +0.5*sigma^2+x(ii)*x(jj));
       %PE contribution
       KPM(ii,jj) = KPM(ii,jj)+(sqrt(pi)/2)*sigma*ev*(((x(ii)-x(jj))/2)^2 ...
                        +0.5*sigma^2);
       
       EM(ii,jj) = EM(ii,jj)+sqrt(pi)*sigma*ev;
   end
end

[eVecs, eVals] = eig(KPM, EM); %Generalized eigenvals
eVals = diag(eVals);
figure; plot(x, (eVecs(:,1:6)).^2);