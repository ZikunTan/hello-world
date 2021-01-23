function [W,E] = bit_forward_msinr_rdrk(K,L,H,V,noise,N_UE,rankU,M,pilot_fwd,pilot_fwd_idv,UE_BS,Ncoop,withNoise,num_reduced_rank,subtract_dirCh)
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
num_reduced_rank = min(num_reduced_rank,N_UE);

for k=1:K
    if rankU(k)==0
        W{k} = zeros(N_UE,1);
    else
        if withNoise
            sigK = H(:,:,k)*HV ...
                + sqrt(noise/2)*randn(N_UE,M)+sqrt(-noise/2)*randn(N_UE,M);   % M number of pilots
            cov_K = sigK*sigK'/M;
            Hv = sigK*pilot_fwd_idv{k}'/M;
        % No noise
        else
            sigK = H(:,:,k)*HV;   % M number of pilots
            cov_K = sigK*sigK'/M+noise*eye(N_UE);
            Hv = sigK*pilot_fwd_idv{k}'/M;
        end
        
        
        %% Reduced rank 
        thisW = [];
        for ns=1:size(Hv,2)
            p = Hv(:,ns);
            
            if subtract_dirCh==1
                cov_ns = cov_K - p*p';
            else
                cov_ns = cov_K;
            end
                
            A = p;
            for n=1:num_reduced_rank-1
                A = [A,cov_ns*A(:,end)];
            end

            cov_rdrk = A'*cov_ns*A;
%             tmpW = A*((cov_rdrk+1e-3*eye(num_reduced_rank))\(A'*p));
            tmpW = A*((cov_rdrk)\(A'*p));
            thisW = [thisW,tmpW];
        end
        
        % Calculate E
        thisE = 1./abs(1-diag(thisW'*Hv));
        for nc=1:size(thisW,2)
            thisW(:,nc) = thisW(:,nc)/norm(thisW(:,nc));
        end
        W{k} = thisW';
        E(k,:) = thisE;
    end
end





