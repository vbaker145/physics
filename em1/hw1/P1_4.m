eps = 0.001;
ri = eps:eps:1-eps;
ro = 1:eps:2;

fo = 1./(ro.^2);
fi1 = ri.*0;
fi2 = ri;
fi31 = 1./ri;
fi32 = ri.^3;

f1 = [fi1 fo];
f2 = [fi2 fo];
f31 = [fi31 fo];
f32 = [fi32 fo];
r = [ri ro];
