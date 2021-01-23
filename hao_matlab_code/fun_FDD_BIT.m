function [rateAll,sinrAll,algs,reswv] = fun_FDD_BIT(simu_param,algs,thisSeed)
% October 2nd, 2018. Hao Zhou
% This is the main function of fdd bit.
% It generates the result for one random drop for algorithms defined in
% 'algs' with same setting/initialization. 
% To average multiple random drops, please run this function multiple times
% with different 'thisSeed'.
[pdb,maxRk,Kp,L,Nr,Nt,S,M,plottestresult,BM_freq_DL,BM_freq_UL,nrank_ue,nrank_bs,rayleigh,rd_tr] = ...
    deal(simu_param.pdb,...
    simu_param.maxRk,...
    simu_param.Kp,...
    simu_param.L,...
    simu_param.Nr,...
    simu_param.Nt,...
    simu_param.S,...
    simu_param.M,...
    simu_param.plottestresult,...
    simu_param.freq_ul_c,...
    simu_param.freq_dl_c,...
    simu_param.nRank_ue,...
    simu_param.nRank_bs,...
    simu_param.rayleigh_fading,...
    simu_param.rd_tr);

subtract_dirCh = simu_param.subtract_dirCh;  % here it is 1
reswv = {};
K = L*Kp;   % total number of user ends
weights = ones(K,1);    % weighted sum rate, K*1 matrix containing ones
% M = max(Nr*K*Ncoop,2000);   % pilot length
% M = 2^11;
Ncoop = 1;  % base station cooperation
% rd_tr = 1;  % random pilots or not


N_UE = Nr;  % wide-band effective antennas  % number of received antennas, 8
N_BS = Nt;                                  % number of transmit antennas, 8
Pmax = 10^((pdb-30)/10);                    % Pmax in unit Watt
Pue = 10^((30-30)/10);                      % unit power each user end ???, in unit Watt
% Pue = Pmax/Kp;

T = 30;    % number of forward-backward iterations
% Nt = 64;    % BS antennas
% Nr = 8;    % UE antennas
% L = 1;      % cells (BSs)
% Kp = 15;  % UE per cell
% pdb = 50; % BS downlink power
% maxRk = 1;  % max rank per UE
randomInit = 1;     % random channel
maxRk_BS = Nt;      % number of transmitting antennas at base station
% Ncoop = 2;    % cooperation cells

noise = 1;
rankU = maxRk*ones(1,K);    % 2 data streams per user, maxRk=2
% S = 12;  % number of paths

%% Initialization users
rng(thisSeed)     % thisSeed=1
% Inter-distance between BSs
minDist = 50;
interBS = 500;
[~,UE_BS,idxBS,AoA,AoD,Dist,Gains,PL] = gen_channel_geo_L(L,Kp,minDist,interBS,Nr,Nt,Ncoop,S,BM_freq_DL,rayleigh);   % rayleigh=1

%% Uplink and downlink bands
F_oneside = 0                                                                                                                                                                          ;

DL_freq = BM_freq_DL+100e3*[-F_oneside:F_oneside];   % BM_freq_DL is downlink central freq

UL_freq = BM_freq_UL+100e3*[-F_oneside:F_oneside];   % BM_freq_UL is uplink central freq

F = length(DL_freq);    % sub-carriers, here it is 1



%% Uplink and downlink channels
H_DL = cell(L,K,F);     % downlink, channel, L number of cells, K total number of user ends, number of sub-carriers
H_UL = cell(L,K,F);     % uplink channel
for f=1:F
    if rayleigh==1                      % if Rayleigh fading applied
        H_rayleigh_tmp = gen_channel_rayleigh(Nt,Nr,PL,K,L);
        H_DL_tmp = H_rayleigh_tmp;
        H_UL_tmp = H_rayleigh_tmp;
    else                                % if Rayleigh fading not applied
        H_DL_tmp = gen_channel_freq(Nt,Nr,S,AoA,AoD,Dist,Gains,DL_freq(f),L,K,BM_freq_DL);
        H_UL_tmp = gen_channel_freq(Nt,Nr,S,AoA,AoD,Dist,Gains,UL_freq(f),L,K,BM_freq_DL);
    end
    
    for l=1:L
        for k=1:K
            H_DL{l,k,f} = H_DL_tmp{l,k};  % downlink channel
            H_UL{l,k,f} = H_UL_tmp{l,k};  % uplink channel
        end
    end
end

% rng('shuffle')
%% Initialize weights
rankBS = zeros(1,L);    % Total initial rank per BS
weightsBS{L} = [];
for k=1:K
    for n=1:Ncoop
        l1 = idxBS(n,k);
        rankBS(l1) = rankBS(l1) + rankU(k);
        weightsBS{l1} = [weightsBS{l1};weights(k)*ones(rankU(k),1)];
    end
end

%% BS V: random initialization

V_FDD = cell(L,F);
for l=1:L
    for f=1:F
%         vini = ones(Nt,rankBS(l));
        vini = sqrt(1/2)*randn(Nt,rankBS(l))+sqrt(-1/2)*randn(Nt,rankBS(l));
        
        vini = vini./vecnorm(vini);
        V_FDD{l,f} = sqrt(Pmax/rankBS(l))*vini;
    end
end

%% combined channel for BIT
H_joint_UE_DL = zeros(N_UE,N_BS*L,K,F);
H_joint_BS_DL = zeros(N_BS,N_UE*K,L,F);
H_joint_UE_UL = zeros(N_UE,N_BS*L,K,F);
H_joint_BS_UL = zeros(N_BS,N_UE*K,L,F);
for f=1:F
    for k=1:K
        HDL = [];
        HUL = [];
        for l=1:L
            HDL = [HDL, H_DL{l,k,f}];
            HUL = [HUL, H_UL{l,k,f}];
        end
        H_joint_UE_DL(:,:,k,f) = HDL;
        H_joint_UE_UL(:,:,k,f) = HUL;
    end
end
for f=1:F
    for l=1:L
        HDL = [];
        HUL = [];
        for k=1:K
            HDL = [HDL, H_DL{l,k,f}'];
            HUL = [HUL, H_UL{l,k,f}'];
        end
        H_joint_BS_DL(:,:,l,f) = HDL;
        H_joint_BS_UL(:,:,l,f) = HUL;
    end
end

%% Generate Pilots
pilot_BS = cell(1,L);          % 1 by L pilots cell
pilot_UE = cell(1,K);          % 1 by K pilots cell
pilot_BS_opt = cell(1,L);
pilot_UE_opt = cell(1,K);
M_tr = max(M,sum(rankU));
% orth_sig = dftmtx(M_tr);
orth_sig = walshcode(M_tr);

ct = 1;
for k=1:K
    rk = rankU(k);
    for n=1:Ncoop
        l = idxBS(n,k);
        if rd_tr==1
            thisPlt = 2*(randi(2,rk,M)-1)-1;
%             thisPlt = sqrt(1/2)*randn(rk,M)+sqrt(-1/2)*randn(rk,M);
            for nrk=1:rk
                thisPlt(nrk,:) = sqrt(M)*thisPlt(nrk,:)/norm(thisPlt(nrk,:));
            end
        else
            thisPlt = orth_sig(ct:ct+rk-1,1:M);
            for nrk=1:rk
                thisPlt(nrk,:) = sqrt(M)*thisPlt(nrk,:)/norm(thisPlt(nrk,:));
            end
        end
        
        thisPlt_opt = orth_sig(ct:ct+rk-1,:);   
        for nrk=1:rk
            thisPlt_opt(nrk,:) = sqrt(M_tr)*thisPlt_opt(nrk,:)/norm(thisPlt_opt(nrk,:));
        end
        
        pilot_BS{l} = [pilot_BS{l};thisPlt];
        pilot_UE{k} = [pilot_UE{k};thisPlt];
        pilot_BS_opt{l} = [pilot_BS_opt{l};thisPlt_opt];
        pilot_UE_opt{k} = [pilot_UE_opt{k};thisPlt_opt];
        ct = ct + rk;
    end
end

disp('pilot_BS:')  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(pilot_BS{1,1})
disp('pilot_UE:')
disp(pilot_UE{1,1})
%% Algorithm comparison
% algs = ["tdd",...
%     "sfdd_try_est_m1_knowDirectCh_minSigProjPower",...
%     "sfdd_try_est_m1_knowDirectCh_angular"];

targetF = ceil(length(DL_freq)/2);  % data sub-carrier
ct = 1;
rateAll = zeros(length(algs),T);
sinrAll = cell(length(algs),1);
while ct<=length(algs)
    rng(thisSeed);
%         algs(ct)
    switch algs(ct)
        
        case {'tdd'}
            % bit solution
            [r,sinr] = fun_fdd_subCases(V_FDD,T,K,L,H_joint_UE_DL,H_joint_UE_UL,H_joint_BS_DL,H_joint_BS_UL,H_DL,noise,...
                N_UE,rankU,idxBS,Ncoop,Pmax,Pue,weights,weightsBS,M,...
                pilot_BS,pilot_UE,UE_BS,S,AoA,AoD,DL_freq,UL_freq,algs(ct),plottestresult,subtract_dirCh);
            
        case {'tdd_opt'}
            % bit solution
            [r,sinr] = fun_fdd_subCases(V_FDD,T,K,L,H_joint_UE_DL,H_joint_UE_UL,H_joint_BS_DL,H_joint_BS_UL,H_DL,noise,...
                N_UE,rankU,idxBS,Ncoop,Pmax,Pue,weights,weightsBS,M_tr,...
                pilot_BS_opt,pilot_UE_opt,UE_BS,S,AoA,AoD,DL_freq,UL_freq,algs(ct),plottestresult,subtract_dirCh);
            
        case {'rdrk'}
            % bit solution
            [r,sinr] = fun_fdd_subCases(V_FDD,T,K,L,H_joint_UE_DL,H_joint_UE_UL,H_joint_BS_DL,H_joint_BS_UL,H_DL,noise,...
                N_UE,rankU,idxBS,Ncoop,Pmax,Pue,weights,weightsBS,M,...
                pilot_BS,pilot_UE,UE_BS,S,AoA,AoD,DL_freq,UL_freq,algs(ct),plottestresult,subtract_dirCh,nrank_ue,nrank_bs);
            
        case {'rdrk_opt'}
            % bit solution
            [r,sinr] = fun_fdd_subCases(V_FDD,T,K,L,H_joint_UE_DL,H_joint_UE_UL,H_joint_BS_DL,H_joint_BS_UL,H_DL,noise,...
                N_UE,rankU,idxBS,Ncoop,Pmax,Pue,weights,weightsBS,M_tr,...
                pilot_BS_opt,pilot_UE_opt,UE_BS,S,AoA,AoD,DL_freq,UL_freq,algs(ct),plottestresult,subtract_dirCh,nrank_ue,nrank_bs);
            
        otherwise
            fprintf('No such algorithm: %s\n',algs(ct))
                        
    end
    
    rateAll(ct,:) = r;
    sinrAll{ct} = sinr;
    ct = ct + 1;
end












