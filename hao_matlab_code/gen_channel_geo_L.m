function [H,UE_perBS,idxBS,AoA,AoD,Dist,Gains,PLs] = gen_channel_geo_L(L,Kp,minDist,interBS,Nr,Nt,Ncoop,S,freq,rayleigh)
% freq = 2.14e9;
d = physconst('LightSpeed')/freq/2;               % half of wavelength
X = [-750,-500,-250,0,250,500,750]'*interBS/500;                  % interBS=500, ' produces complex conjugate transpose
Y = [0,250*sqrt(3),0,250*sqrt(3),0,250*sqrt(3),0]'*interBS/500;

BS_loc = [X(1:L),Y(1:L)];     % L=7
itf_region = BS_loc([1:L,1],:);
lb = min(min(BS_loc))-interBS/2;  % -750-250
ub = max(max(BS_loc))+interBS/2;  % 750+250
% BS_loc = [0,0];

K = Kp*L;               % total number of user ends
UE_perBS = cell(1,L);   % 1 by L cell array of empty matrices
UE_pos = [];
idxBS = ones(Ncoop,K);  % Ncoop by K matrix, 1 by 70 matrix
kk = 1;

for l=1:L               % for each cell, L=7
    ct = 1;
    while ct<=Kp        % for users in each cell, Kp=10
        uepos = (ub-lb)*rand(1,2)+lb;    % 1 by 2 random matrix
%         if L>2
%             if ~inpolygon(uepos(1),uepos(2),itf_region(:,1),itf_region(:,2))
%                 continue;
%             end
%         end
        [md,mdid] = min(vecnorm((uepos-BS_loc)'));
        if mdid~=l || md<=minDist || md>interBS/2
            continue;
        end
        idxBS(kk) = l;
        UE_pos = [UE_pos; uepos];
        UE_perBS{l} = [UE_perBS{l},kk];
        ct = ct + 1;
        kk = kk + 1;
    end
end

% circle(0,0,1)
% plot(BS_loc(:,1),BS_loc(:,2),'x')
% hold on;
% plot(UE_pos(:,1),UE_pos(:,2),'o')

AoA = cell(L,K);
AoD = cell(L,K);
Dist = cell(L,K);
Gains = cell(L,K);
H = cell(L,K);
PLs = zeros(L,K);
for k=1:K
    for l=1:L
        thidD = norm(BS_loc(l,:)-UE_pos(k,:));    % large-scale fading
        thisPL = 32.4+20*log10(3.5)+31.9*log10(thidD) - 174 + 10*log10(1e6) +8*randn;
        thisPL_mag = sqrt(10^(-thisPL/10));
        PLs(l,k) = thisPL_mag;
        if rayleigh==0
            [H{l,k},AoA{l,k},AoD{l,k},Dist{l,k},Gains{l,k}] = fun_geo_ch(Nt,Nr,S,BS_loc(l,:),UE_pos(k,:),freq,d,thisPL_mag);
        end
    end
end























