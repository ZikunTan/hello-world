function H = gen_channel_rayleigh(Nt,Nr,PL,K,L)
H = cell(L,K);  % number of cell by number of total user ends, L by K

for k=1:K
    for l=1:L
        H{l,k} = sqrt(PL(l,k)/2)*randn(Nr,Nt)+sqrt(-PL(l,k)/2)*randn(Nr,Nt);   % generate Rayleigh fading channel                                                                          % PL(l,k)
    end                                                                        % PL(l,k) is the mean power of channel impulse response
end
