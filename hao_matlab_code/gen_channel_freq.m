function H = gen_channel_freq(N,M,S,AoA,AoD,Dist,Gains,freq,L,K,BM_freq)
H = cell(L,K);
lambda = physconst('LightSpeed')/freq;
lambda0 = physconst('LightSpeed')/freq;

d = physconst('LightSpeed')/BM_freq/2;

for k=1:K
    for l=1:L
        thisH = zeros(M,N);
        for s=1:S
            thisH = thisH + Gains{l,k}(s)*exp(-1j*2*pi*Dist{l,k}(s)/lambda0)...
                *exp(-1j*2*pi*d/lambda*cos(AoA{l,k}(s))*[0:M-1]')...
                *exp(1j*2*pi*d/lambda*cos(AoD{l,k}(s))*[0:N-1]);
        end
        H{l,k} = thisH;
    end
end

