clear all; close all

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
        if (l==lp)
           if (mp==(m+1))
               %L+
               hlx(ridx, cidx) = hlx(ridx, cidx)+sqrt( l*(l+1)-m*(m+1) );
           end
           if (mp==(m-1))
               %L- 
               hlx(ridx, cidx) = hlx(ridx, cidx)+sqrt( l*(l+1)-m*(m-1) );
           end
        end
        
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
hlx = hlx./2; %Make Lx from L+, L-

%Use eig to find eigenvalues/vectors

ev0 = eig(h0);
evz = eig(hlz);
evx = eig(hlx);
eve = eig(he);

figure; plot(evx,'go')
hold on; plot(evx, 'rx')
title('Spectrum of L_z and L_x')
xlabel('Eigenvalue #'); ylabel('Eigenvalue')

figure; plot(eve, 'x')
title('Spectrum of H_{E}')
he = he./emFac; %Make Stark effect same magitude as Zeeman

%Find total system spectrum
ndegs = 100;
theta = linspace(0, pi, ndegs);
for didx = 1:ndegs
    ht = h0+he+hlz*cos(theta(didx))+hlx*sin(theta(didx));
    evt(:,didx)=eig(ht);
end

figure; plot(theta.*180/pi, evt'); 
title('Energy spectrum of sodium valence electron');
xlabel('theta (degrees)'); ylabel('energy')