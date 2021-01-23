function [W,E] = bit_forward_msinr(K,L,H,V,noise,N_UE,rankU,M,pilot_fwd,pilot_fwd_idv,UE_BS,Ncoop,withNoise,subtract_dirCh)
% This function does the forward training 
% Output the UE's combiners W

N_BS = size(V{1},1);
W = cell(1,K);
maxRank = max(rankU);
E = ones(K,maxRank*Ncoop);
HV = zeros(N_BS,M,L);

for l=1:L
    if ~isempty(UE_BS{l})
        HV(:,:,l) = V{l}*pilot_fwd{l};
    else
        HV(:,:,l) = zeros(N_BS,M);
    end
end

HV = reshape(permute(HV,[1,3,2]),N_BS*L,M);
for k=1:K
    if rankU(k)==0
        W{k} = zeros(N_UE,1);
    else
        if withNoise
            sigK = H(:,:,k)*HV ...
                + sqrt(noise/2)*randn(N_UE,M)+sqrt(-noise/2)*randn(N_UE,M);   % M number of pilots
            cov_K = sigK*sigK'/M;
            Hv = sigK*pilot_fwd_idv{k}'/M;
        else
        % No noise
            sigK = H(:,:,k)*HV;   % M number of pilots
            cov_K = sigK*sigK'/M+noise*eye(N_UE);
            Hv = sigK*pilot_fwd_idv{k}'/M;
        end
        
        if subtract_dirCh==1
            thisW = [];
            for ns=1:size(Hv,2)
                thisW = [thisW,(cov_K-Hv(:,ns)*Hv(:,ns)')\Hv(:,ns)];
            end
        else
            thisW = cov_K\Hv;
        end
        % Calculate E
        thisE = 1./abs(1-diag(thisW'*Hv));
        for ns=1:size(thisW,2)
            thisW(:,ns) = thisW(:,ns)/norm(thisW(:,ns));
        end
        W{k} = thisW';
        E(k,:) = thisE;
    end
end





