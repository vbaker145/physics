clear all; close all

%Set up mass/flacor mixing matrix
t12 = 0.5*(asin(sqrt(0.861)));
t13 = 0.5*(asin(sqrt(0.092)));
t23 = 0.5*(asin(sqrt(0.97)));

s12 = sin(t12); c12 = cos(t12);
s13 = sin(t13); c13 = cos(t13);
s23 = sin(t23); c23 = cos(t23);

ut1 = [1 0 0; 0 c23 s23; 0 -s23 c23];
ut2 = [c13 0 s13; 0 1 0; -s13 0 c13];
ut3 = [c12 s12 0; -s12 c12 0; 0 0 1];

ut = ut1*ut2*ut3;

%Find inverse to resolve flavor (1) state to mass states
%(Apparently I don't beleive this transformation is unitary, so instead of
%just transposing it I'll let MATLAB calculate the inverse)
uti = inv(ut);

%Calculate the mass state mixture at t=0
m_f1 = uti*[1; 0; 0];

%Propagate the energy states forward in time from t=0
t=linspace(0,4*pi, 3000);
m1 = [1 0 0]*m_f1*exp(-1i*t);
m2 = [0 1 0]*m_f1*exp(-4i*t);
m3 = [0 0 1]*m_f1*exp(-24i*t);

%Calculate flavor mixtures from mass states
f = ut*[m1; m2; m3];
f1 = f(1,:);
f2 = f(2,:)
f3 = f(3,:);

%Calculate flavor probabilities
p100 = f1.*conj(f1);
p010 = f2.*conj(f2);
p001 = f3.*conj(f3);

%Plot amplitude
figure; plot(real(f1))
hold on; plot(real(f2), 'r')
hold on; plot(real(f3),'g')

%Plot probability
figure; plot(p100)
hold on; plot(p010, 'r')
hold on; plot(p001, 'g')
xlabel('t'); ylabel('probability')
legend('P(flavor1)', 'P(flavor2)', 'P(flavor3)');