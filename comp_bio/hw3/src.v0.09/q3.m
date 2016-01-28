%f(x,y) = x^2+2y^2
[x y] = meshgrid(-20:0.1:20, -20:0.1:20);
figure(10); surf(x.^2+2.*(y.^2));

x0 = 5; y0 = 4;

eps = 1e-6;
err = x0^2+2*y0^2;
pos = [x0 y0 err];

while err > eps
   %Calculate direction
   dir_norm = sqrt(4*x0^2+16*y0^2);
   dir = [-2*x0 -4*y0 ]./dir_norm;
   
   %Find slope-intercept of search line
   m = dir(2)/dir(1);
   intercept = y0-x0*m;
   
   %Calculate polynomial coefficients, find derivative equation = 0
   a = 2*(1+2*m^2);
   b = 2*2*m*intercept;
  
   %New x0/y0
   x0 = -b/a;
   y0 = m*x0+intercept;
   
   err = x0^2+2*y0^2; 
   pos = [pos; x0 y0 err];
end