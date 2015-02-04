
x=linspace(-10,10,100);
a = [.01,.1,1,10, 10e6];
figure(10)
for jj=1:length(a)
    c1 = x+sign(x).*(x.^2/(2*a(jj)));
    hold on; plot(x, c1)
    hold on; plot(x,-c1, 'r')
end