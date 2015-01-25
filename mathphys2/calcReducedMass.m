function [ rm ] = calcReducedMass(mep, mpp)
%Reduced mass for hydrogen-like systems
%   mep - mass of electron-like particle, kg
%   mpp - mass of proton-like particle, kg

me = 9.11e-31;
mp = 1.67e-27;
ea = mep/me; %Scale factor for electron particle
pa = mpp/mp; %Scale factor for proton particle

rm = ea*pa*((me+mp)/(ea*me+pa*mp));

end

