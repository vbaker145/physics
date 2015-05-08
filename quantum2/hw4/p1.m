t=0:pi/100:2*pi;

s1 = 5/9+(2/9)*(cos(6*t));

s2 = 4/9+(2/9)*(cos(6*t));

ent = -s1.*log(s1)-s2.*log(s2);

figure; plot(t./(2*pi), ent)

