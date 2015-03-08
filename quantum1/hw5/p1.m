clear all

N = 3;
nc = 9;
emFac = 4.5; %Porpotionality factor between Stark/Zeeman effect energy levels

h0  = zeros(nc,nc);
hlz = zeros(nc,nc);
hlx = zeros(nc,nc);
he  = zeros(nc,nc);

lvals = [0 1 1 1 2 2 2 2 2];
mvals = [0 1 0 -1 2 1 0 -1 -2];
lidx = 0;
for cidx = 1:nc
    for ridx = 1:nc
        l = lvals(cidx);  m = mvals(cidx);
        lp = lvals(ridx); mp = mvals(ridx);
        
        %H0
        if cidx == ridx
           if cidx == 1
               h0(ridx, cidx) = h0(ridx, cidx)+0;
           elseif (1<cidx) &&(cidx< 5)
               h0(ridx, cidx) = h0(ridx, cidx)+2;
           else
               h0(ridx, cidx) = h0(ridx, cidx)+3;
           end
        end
        
        %Lz
        if (l==lp) && (m==mp) 
           hlz(ridx, cidx) = hlz(ridx,cidx) + m; 
        end
        %Lx
        
        %E.r
        if (l-lp) == 1
           if (mp == m)
              angFac = sqrt( ((l+m)*(l-m)) / ((2*l+1)*(2*l-1)) );
              rFac   = (3/2)*N*sqrt(N^2-l^2);
              he(ridx,cidx) = he(ridx,cidx)+rFac*angFac;
           end
        end
        if (lp-l) == 1
           if (mp == m)
              angFac = sqrt( ((lp+mp)*(lp-mp)) / ((2*lp+1)*(2*lp-1)) );
              rFac   = (3/2)*N*sqrt(N^2-lp^2);
              he(ridx,cidx) = he(ridx,cidx)+rFac*angFac;
           end
        end
        
    end
end

he = he./emFac; %Make Stark effect same magitude as Zeeman

ht = h0+hlz+hlx+he;

%Use eig to find eigenvalues/vectors