clear all;                         % clear all variables, functions and MEX files
close all;                         % close all figure windows
clc                                % clear all content in command window
warning off;                       % not display warning information

simu_param.pdb = 35;               % power, dbm
simu_param.Kp = 10;                % number of users per cell
simu_param.L = 7;                  % number of cells
simu_param.Nr = 8;                 % number of antennas at mobile
simu_param.Nt = 8;                 % number of antennas at base station
simu_param.S = 24;                 % number of paths
simu_param.freq_ul_c = 680e6;      % uplink central frequency
simu_param.freq_dl_c = 680e6;      % downlink central frequency
simu_param.plottestresult = 0;     % plot sum rate vs F-B iteration
simu_param.rayleigh_fading = 1;
simu_param.rd_tr = 0;              % random training pilots or orthogonal (binary)
simu_param.subtract_dirCh = 1;     % subtract direct channel 

%% Algorithms
%   1. Full rank filter w/ training (tdd): max-sinr w/ equal power
%   2. Full rank filter w/ optimization (tdd_opt): F-B optimization
%   3. Reduced rank w/ training (rdrk): max-sinr w/ equal power
%   4. Reduced rank w/ optimization (rdrk_opt): F-B optimization w/ reduced rank
algs = ["tdd","tdd_opt","rdrk","rdrk_opt"];

thisSeed = 1;           % random drop seed

%% Default parameters
simu_param.maxRk = 2;              % number of data streams per user
simu_param.M = 512;                % number of pilots
simu_param.nRank_ue = 4;           % number of reduced rank at user end
simu_param.nRank_bs = 4;           % number of reduced rank at base station

% Save default parameters to param
param_default = simu_param;

%% Example 1: Plot sum rate w.r.t. rank of reduced rank filter
[simu_param.M,simu_param.nRank_ue,simu_param.nRank_bs,simu_param.maxRk] = deal(param_default.M,param_default.nRank_ue,param_default.nRank_bs,param_default.maxRk);    % assign the right side to the left side
% Initialization
nRank_all = 1:8;       % filter rank to be tested, from 1 to 8, increment 1
[res,sinrAll,~] = fun_FDD_BIT(simu_param,algs,thisSeed);    % only need variables res and sinrALL
res_all = repmat(res(:,end),1,length(nRank_all));
% Run for different test cases
for n=1:length(nRank_all)
    simu_param.nRank_ue = nRank_all(n);
    simu_param.nRank_bs = nRank_all(n);
    fprintf('Reduced rank ue/bs: %d,%d\n',simu_param.nRank_ue,simu_param.nRank_bs);
    [res,sinrAll,~] = fun_FDD_BIT(simu_param,["rdrk","rdrk_opt"],thisSeed);
    % res: size N_algs x F-B iterations
    res_all(3:4,n) = res(:,end);
end
% Plot results
figure;
plot(nRank_all,res_all','-o');
legend(algs);
xlabel('Reduced rank')
ylabel('Sum rate, bps/Hz')
title(sprintf('L%dNt%dNr%dP%dK%d,nstream%dM%d',simu_param.L,simu_param.Nt,simu_param.Nr,simu_param.pdb,simu_param.Kp,simu_param.maxRk,simu_param.M))

%% Example 2: Plot sum rate w.r.t. training length
[simu_param.M,simu_param.nRank_ue,simu_param.nRank_bs,simu_param.maxRk] = deal(param_default.M,param_default.nRank_ue,param_default.nRank_bs,param_default.maxRk);
% Initialization
M_all = 2.^(4:10);
[res,sinrAll,~] = fun_FDD_BIT(simu_param,algs,thisSeed);
res_all = repmat(res(:,end),1,length(M_all));
% Run for different test cases
for n=1:length(M_all)
    simu_param.M = M_all(n);
    fprintf('Training length: %d\n',simu_param.M);
    [res,sinrAll,~] = fun_FDD_BIT(simu_param,["tdd","rdrk"],thisSeed);
    % res: size N_algs x F-B iterations
    res_all([1,3],n) = res(:,end);
end
% Plot results
figure;
plot(M_all,res_all','-o');
legend(algs);
xlabel('Training length')
ylabel('Sum rate, bps/Hz')
title(sprintf('L%dNt%dNr%dP%dK%d,nRank%dnstream%d',simu_param.L,simu_param.Nt,simu_param.Nr,simu_param.pdb,simu_param.Kp,simu_param.nRank_bs,simu_param.maxRk))

%% Example 3: Plot sum rate w.r.t. precoder rank
[simu_param.M,simu_param.nRank_ue,simu_param.nRank_bs,simu_param.maxRk] = deal(param_default.M,param_default.nRank_ue,param_default.nRank_bs,param_default.maxRk);
% Initialization
maxRk_all = 1:4;
res_all = zeros(4,length(maxRk_all));
% Run for different test cases
for n=1:length(maxRk_all)
    simu_param.maxRk = maxRk_all(n);
    fprintf('Precoder rank: %d\n',simu_param.maxRk);
    [res,sinrAll,~] = fun_FDD_BIT(simu_param,algs,thisSeed);
    % res: size N_algs x F-B iterations
    res_all(:,n) = res(:,end);
end
% Plot results
figure;
plot(maxRk_all,res_all','-o');
legend(algs);
xlabel('Precoder rank')
ylabel('Sum rate, bps/Hz')
title(sprintf('L%dNt%dNr%dP%dK%d,nRank%dM%d',simu_param.L,simu_param.Nt,simu_param.Nr,simu_param.pdb,simu_param.Kp,simu_param.nRank_bs,simu_param.M))

            


