function rate = calcSinr(W,V,H,K,L,rankU,noise,idxBS,Ncoop)
% This function calculates the downlink sinr

rate = zeros(K,max(rankU)*Ncoop);
ct = ones(1,L);
for k=1:K
    % Received signals
    sig_all = [];
    for li=1:L
        sig_all = [sig_all, H{li,k}*V{li}];
    end
    rk = rankU(k);
    for n=1:Ncoop
        l = idxBS(n,k);
        for r=1:rk
            numr = norm(W{k}((n-1)*rk+r,:)*H{l,k}*V{l}(:,ct(l)))^2;
            ct(l) = ct(l) + 1;
            d_tmp = W{k}((n-1)*rk+r,:)*sig_all;
            denm = sum(abs(d_tmp).^2)-numr+noise*norm(W{k}((n-1)*rk+r,:))^2;
            this_sinr = numr/denm;
            if isnan(this_sinr)
                rate(k,(n-1)*rk+r) = 0;
            else
                rate(k,(n-1)*rk+r) = this_sinr;
            end
        end
    end
end

