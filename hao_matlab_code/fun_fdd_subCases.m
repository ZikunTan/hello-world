function [rate,sinr] = fun_fdd_subCases...
    (V_bit,T,K,L,H_joint_UE_DL,H_joint_UE_UL,H_joint_BS_DL,H_joint_BS_UL,H_DL,...
    noise,N_UE,rankU,idxBS,Ncoop,Pmax,Pue,weights,weightsBS,M,pilot_BS,pilot_UE,...
    UE_BS,S,AoA,AoD,DL_freq,UL_freq,algType,plottestresult,subtract_dirCh,varargin)

% October 3rd, 2018. Hao Zhou
% This function runs sub-cases of fdd bit.
% Supporting cases:
%   1. TDD (tdd): max-sinr w/ equal power
%   2. tdd_opt: no training noise (F-B optimization)
%   3. Reduced rank (rdrk): max-sinr w/ equal power
%   4. rdrk_opt: no training noise (F-B optimization w/ reduced rank)

withNoise = 1;
if ~isempty(varargin)
    nrank_ue = varargin{1};
    nrank_bs = varargin{2};    
end

if plottestresult
    disp(algType);
end
r_bit = zeros(1,T); % rate over iterations
F = length(DL_freq);
targetF = ceil(length(DL_freq)/2);  % data sub-carrier
estAoA = cell(1,T);
estAoD = cell(1,T);

res = cell(2,T);
sinr_overFreq = zeros(K,F,T);
for t=1:T
    if plottestresult
        fprintf('--- F-B iteration: %d\n',t);
    end
    preV = V_bit;
    switch algType
        case 'tdd'        
            [W_bit,~] = bit_forward_msinr(K,L,H_joint_UE_DL(:,:,:,targetF),V_bit,noise,N_UE,rankU,M,pilot_BS,pilot_UE,UE_BS,Ncoop,withNoise,subtract_dirCh);
            V_bit = bit_backward_msinr(K,L,N_UE,rankU,W_bit,H_joint_BS_DL(:,:,:,targetF),Pue,Pmax,V_bit,noise,M,pilot_BS,pilot_UE,weights,weightsBS,Ncoop,withNoise,subtract_dirCh);
           
        case 'tdd_opt'          
            withNoise = 0;
            [W_bit,~] = bit_forward_msinr(K,L,H_joint_UE_DL(:,:,:,targetF),V_bit,noise,N_UE,rankU,M,pilot_BS,pilot_UE,UE_BS,Ncoop,withNoise,subtract_dirCh);
            V_bit = bit_backward_msinr(K,L,N_UE,rankU,W_bit,H_joint_BS_DL(:,:,:,targetF),Pue,Pmax,V_bit,noise,M,pilot_BS,pilot_UE,weights,weightsBS,Ncoop,withNoise,subtract_dirCh);
            
        case 'rdrk'
            [W_bit,~] = bit_forward_msinr_rdrk(K,L,H_joint_UE_DL(:,:,:,targetF),V_bit,noise,N_UE,rankU,M,pilot_BS,pilot_UE,UE_BS,Ncoop,withNoise,nrank_ue,subtract_dirCh);
            V_bit = bit_backward_msinr_rdrk(K,L,N_UE,rankU,W_bit,H_joint_BS_DL(:,:,:,targetF),Pue,Pmax,V_bit,noise,M,pilot_BS,pilot_UE,weights,weightsBS,Ncoop,withNoise,nrank_bs,subtract_dirCh);
            
        case 'rdrk_opt'
            withNoise = 0;
            [W_bit,~] = bit_forward_msinr_rdrk(K,L,H_joint_UE_DL(:,:,:,targetF),V_bit,noise,N_UE,rankU,M,pilot_BS,pilot_UE,UE_BS,Ncoop,withNoise,nrank_ue,subtract_dirCh);
            V_bit = bit_backward_msinr_rdrk(K,L,N_UE,rankU,W_bit,H_joint_BS_DL(:,:,:,targetF),Pue,Pmax,V_bit,noise,M,pilot_BS,pilot_UE,weights,weightsBS,Ncoop,withNoise,nrank_bs,subtract_dirCh);
            
   
        otherwise
            disp('No such algorithm')
    end
    
        
    if size(preV,2)>1
        tmpV = cell(1,L);
        for l=1:L
            tmpV{l} = preV{l,targetF};
        end
        sinr_jt = calcSinr(W_bit,tmpV,H_DL(:,:,targetF),K,L,rankU,noise,idxBS,Ncoop);
    else
        sinr_jt = calcSinr(W_bit,preV,H_DL(:,:,targetF),K,L,rankU,noise,idxBS,Ncoop);
    end
    r_bit(t) = sum(sum(log2(1+sinr_jt)));
    

    if strcmp(algType,'rbf')
        r_bit(t:end) = r_bit(t);
        break;
    end
    if t>10 && abs((r_bit(t)-r_bit(t-1))/r_bit(t-1))<1e-3
        r_bit(t:end) = r_bit(t);
        break;
    end
    
end
if plottestresult
    hold on;plot(r_bit);
end
rate = r_bit;
sinr = 10*log10(sinr_jt);

