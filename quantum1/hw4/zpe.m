function [ zpe ] = zpe( N )
%Calculate zero-point energy of N coupled harmonic oscillators
theta = (pi/2)/(N+1);
zpe = (cos(theta)+sin(theta)-1)/(1-cos(theta));

end

