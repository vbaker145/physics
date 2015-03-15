clear all; close all

%Set up mass/flavor mixing matrix
theta = asin(.1);
st = sin(theta); ct = cos(theta);
ut = [ct st; -st ct];

%Find inverse to resolve flavor (1) state to mass states
%(Apparently I don't beleive this transformation is unitary, so instead of
%just transposing it I'll let MATLAB calculate the inverse)
uti = inv(ut);

%Find mass mixture for flavor [1 0]
m_f1 = uti*[1; 0];

%Propagate mass (energy) states forward in time
t=linspace(0,10*pi, 3000);
m1 = [1 0]*m_f1*exp(-1i*t);
m2 = [0 1]*m_f1*exp(-10i*t);

%Find flavor mixtures from mass mixtures
f = ut*[m1; m2];
f1 = f(1,:);
f2 = f(2,:)
p10 = f1.*conj(f1);
p01 = f2.*conj(f2);

%figure; plot(abs(f1))
%hold on; plot(abs(f2), 'r')

figure; plot(p10)
hold on; plot(p01, 'r')

