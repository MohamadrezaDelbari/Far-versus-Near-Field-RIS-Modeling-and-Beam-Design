%%----------------------------------------------------%%
%%----- Mohamadreza Delbari
%%      Please cite our paper:
%%----- DOI: https://arxiv.org/pdf/2401.08237
%%----------------------------------------------------%%
function [W] = near_ana_2D(f,Ny,Nz,pr,pt,pirs,ppirs)
c=3*10^8;
lamda=c/f;
dy=lamda/2;
dz=lamda/2;
k=2*pi/lamda;
N=Ny*Nz;
n_z = Nz-1:-1:0;
n_y = Ny-1:-1:0;

%______________________Position________________________
prx=linspace(pr(1,1),pr(2,1),N); %Receiver
pry=linspace(pr(1,2),pr(2,2),N);


dt=vecnorm([pt(1,1) pt(1,2) pt(1,3)]'-ppirs');
dr=vecnorm([prx' pry' pr(1,3)*ones(N,1)]'-ppirs');
dn=vecnorm(pirs'-[prx' pry' pr(1,3)*ones(N,1)]');

w = exp(1i*k*(dr-dn-dt));
W=diag(w);
%save("near_ana_1D.mat","w","W")
end