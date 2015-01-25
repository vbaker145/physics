function hp = genhermite(n,x,b)
if n==1
    hp = 1
    return
end
if n==2
    hp = 2*x*b
    return
end

hp = 2*b*x;
hpm1 = 1;
for jj=1:n-2
    tmp = hp;
    hp = 2*b*x*hp-2*jj*hpm1;
    hpm1 = tmp;
end

hp = simplify(hp);


