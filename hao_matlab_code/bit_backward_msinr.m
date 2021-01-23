function V = bit_backward_msinr(K,L,N_UE,rankU,W,H,Pue,Pmax,pre_V,noise,M,pilot_BS,pilot_UE,alpha,alphaBS,Ncoop,withNoise,subtract_dirCh)
for k=1:K
    W{k} = W{k}*sqrt(Pue/size(W{k},1));
end
N_BS = size(H(:,:,1),1);
V = cell(L,1);
HV = zeros(N_UE,M,K);
for k=1:K
    if rankU(k)>0
        HV(:,:,k) = sqrt(alpha(k))*W{k}'*pilot_UE{k};
    end
end
HV = reshape(permute(HV,[1,3,2]),N_UE*K,M);

for l=1:L
    NrankBS = size(pre_V{l},2);
    if NrankBS>0
        pV = pre_V{l};
        % Uplink received signal
        if withNoise
            sigL = H(:,:,l)*HV ...
                + sqrt(noise/2)*randn(N_BS,M)+sqrt(-noise/2)*randn(N_BS,M);   % M number of pilots
            cov_A = sigL*sigL'/M;    % Signal covariance matrix
            Hv = sigL*pilot_BS{l}'/M;
        % No noise
        else
            sigL = H(:,:,l)*HV;   % M number of pilots
            cov_A = sigL*sigL'/M+noise*eye(N_BS);    % Signal covariance matrix
            Hv = sigL*pilot_BS{l}'/M;
        end
        
        if subtract_dirCh==1
            thisV = [];
            for ns=1:size(Hv,2)
                thisV = [thisV,(cov_A-Hv(:,ns)*Hv(:,ns)')\Hv(:,ns)];
            end
        else
            thisV = cov_A\Hv;
        end
        
        for nc=1:size(thisV,2)
            thisV(:,nc) = thisV(:,nc)/norm(thisV(:,nc));
        end
        V{l} = thisV*sqrt(Pmax/size(thisV,2));
    else
        V{l} = zeros(N_BS,0);
    end
end
