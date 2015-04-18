clear all

jx2 = (1/4)*[2 0 2; 0 4 0; 2 0 2 ];
jy2 = (1/4)*[2 0 -2; 0 4 0; -2 0 2];
jz2 = [1 0 0; 0 0 0; 0 0 1];

%Commutation relations
cxy = jx2*jy2-jy2*jx2;
cyz = jy2*jz2-jz2*jy2;
czx = jz2*jx2-jx2*jz2;

%Sum of all three squares
s = jx2+jy2+jz2;
