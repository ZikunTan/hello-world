function [H,AoA,AoD,Dist,gains] = fun_geo_ch(N,M,S,posBS,posUE,freq,d,thisPL_mag)
H = zeros(M,N);
AoD = atan2((posUE(2)-posBS(2)),(posUE(1)-posBS(1)));
candidateangles = acos(-1:0.0002:1-0.02);
% candidateangles = acos(-0.9:0.02:0.9-0.02);
% if posUE(2)>0
%     if AoD<0
%         AoD = pi+AoD;
%     end
% else
%     if AoD>0
%         AoD = AoD-pi;
%     end
% end 
AoA = mod(pi+AoD,2*pi);
Dist = norm(posUE-posBS);
for s=2:S
    Dist = [Dist, norm(posUE-posBS) + rand*10];
end

% AoA(2:S) = acos(mod(linspace(0,2-2/S,S-1)*(rand*0.5+0.5),2)-1);
% AoD(2:S) = acos(mod(linspace(0,2-2/S,S-1)*(rand*0.5+0.5),2)-1);

[~,id] = min(abs(cos(AoA)-cos(candidateangles)));
Ndiff = floor(length(candidateangles)/S);
aoa_ids = [id,id+(1:S-1).*Ndiff+randi( Ndiff,1,S-1)];
[~,id] = min(abs(cos(AoD)-cos(candidateangles)));
aod_ids = [id,id+10+(0:S-2).*Ndiff+randi( Ndiff,1,S-1)];

AoA = candidateangles(mod(aoa_ids,length(candidateangles))+1);
AoD = candidateangles(mod(aod_ids,length(candidateangles))+1);
for s=1:S
    if abs(cos(AoA(s))+1)<1e-3
        AoA(s) = acos(-1+1e-3);
    end
    if abs(cos(AoA(s))-1)<1e-3
        AoA(s) = acos(1-1e-3);
    end
    if abs(cos(AoD(s))+1)<1e-3
        AoD(s) = acos(-1+1e-3);
    end
    if abs(cos(AoD(s))-1)<1e-3
        AoD(s) = acos(1-1e-3);
    end
end

lambda = physconst('LightSpeed')/freq;
gains = randn(1,S) + 1j*randn(1,S);
gains = sort(abs(gains),'descend');
% gains(2:end) = gains(2:end)/2;
gains = thisPL_mag*gains/gains(1);
for s=1:S
    H = H + gains(s)*exp(-1j*2*pi*Dist(s)/lambda)...
        *exp(-1j*2*pi*d/lambda*cos(AoA(s))*[0:M-1]')...
        .*exp(1j*2*pi*d/lambda*cos(AoD(s))*[0:N-1]);
end

