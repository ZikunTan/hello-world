X = [-750,-500,-250,0,250,500,750]'*500/500;                  % interBS=500, ' produces complex conjugate transpose
Y = [0,250*sqrt(3),0,250*sqrt(3),0,250*sqrt(3),0]'*500/500;
L=6
BS_loc = [X(1:L),Y(1:L)];
display(BS_loc)
itf_region = BS_loc([1:L,1],:)
display(itf_region)

min(min(BS_loc))