%Calculate binding energy for hydrogen-like systems
m_e = 9.11e-31;
m_p = 1.67e-27;
be_h = -13.6;
d_h  = 1.058;

mev2kg = 1.78e-30;
m_muon = 105.66*mev2kg;
m_pion = 139.57*mev2kg;

%Hydrogen (test)
be_h2 = calcReducedMass(m_e, m_p)*be_h
d_h2  = (1/calcReducedMass(m_e, m_p))*d_h

%Helium (-1 electron)
be_hel = calcReducedMass(m_e, 2*m_p)*be_h
d_hel  = (1/calcReducedMass(m_e, 2*m_p))*d_h

%Muon/pion sub for electron
be_mu = calcReducedMass(m_muon, m_p)*be_h
d_mu  = (1/calcReducedMass(m_muon, m_p))*d_h
be_pi = calcReducedMass(m_pion, m_p)*be_h
d_pi  = (1/calcReducedMass(m_pion, m_p))*d_h

%Positronium
be_ee = calcReducedMass(m_e, m_e)*be_h
d_ee  = (1/calcReducedMass(m_e, m_e))*d_h

%Muonium and positronium
be_muonium = calcReducedMass(m_muon, m_muon)*be_h
d_muonium  = (1/calcReducedMass(m_muon, m_muon))*d_h
be_pionium = calcReducedMass(m_pion, m_pion)*be_h
d_pionium  = (1/calcReducedMass(m_pion, m_pion))*d_h

%Silicon exciton
si_m = calcReducedMass(0.8*m_e, 0.4*m_e);
be_si = be_h*si_m*1/(11.9^2);
d_si = d_h*(1/si_m)*11.9;

%GaAs exciton
gaas_m = calcReducedMass(0.07*m_e, 0.4*m_e);
be_gaas = be_h*gaas_m*(1/12.5^2);
d_gaas = d_h*(1/gaas_m)*12.5;

