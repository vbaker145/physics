clear all;

vals = [1,0,-1];
Nvals = 3;
N=3;

state = [];
for kk=1:length(vals)
    %All are first val
    idx = kk;
    state = [state; vals(kk) vals(kk) vals(kk)];
    
    %Two first val
    idx= incMod(kk);
    state = [state; vals(idx) vals(kk) vals(kk); vals(kk) vals(idx) vals(kk); vals(kk) vals(kk) vals(idx)];
    idx= incMod(idx);
    state = [state; vals(idx) vals(kk) vals(kk); vals(kk) vals(idx) vals(kk); vals(kk) vals(kk) vals(idx)];

    %One first val, other two second or third
    idx= incMod(kk);
    state = [state; vals(kk) vals(idx) vals(idx); vals(idx) vals(kk) vals(idx); vals(idx) vals(idx) vals(kk) ];
    idx= incMod(idx);
    state = [state; vals(kk) vals(idx) vals(idx); vals(idx) vals(kk) vals(idx); vals(idx) vals(idx) vals(kk) ];
    
    %one first val, one each from other two
    state = [state; vals(kk) vals(incMod(kk)) vals(incMod(incMod(kk)))];
    state = [state; vals(kk) vals(incMod(incMod(kk))) vals(incMod(kk))];
end

state = unique(state, 'rows');

sv = sum(state');

figure; hist(sv, [-3:3]);

%Calculate entropy, all states same
Sall=0; Ns= length(sv);
for jj=1:Ns
   Sall = Sall - (1/Ns)*log(1/Ns);
end

%Calculate entropy, only M=0
S0=0; N0 = length(find(sv==0));
for jj=1:N0
   S0 = S0 - (1/N0)*log(1/N0);
end
%Calculate entropy, only M=1
S1=0; N1 = length(find(sv==1));
for jj=1:N1
   S1 = S1 - (1/N1)*log(1/N1);
end

%Calculate entropy, only M=3
S3=0; N3 = length(find(sv==3));
for jj=1:N3
   S3 = S3 - (1/N3)*log(1/N3);
end

