
clear all; close all

va = 9e-6; na = 3;
vb = 4e-6; nb = 2;
et = []
fracs = 0:0.01:1;
for eaf = fracs
    ea = eaf*80;
    eb = (1-eaf)*80;
    et =  [et (na*va*ea)^(1/3)+(nb*vb*eb)^(1/3)];
end
figure; plot(fracs, et)

ebea = sqrt((nb*vb)/(na*va))