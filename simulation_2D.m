%%----------------------------------------------------%%
%%----- Mohamadreza Delbari
%%      Please cite our paper:
%%----- DOI: https://arxiv.org/pdf/2401.08237
%%----------------------------------------------------%%

% Delete all the path cashes
restoredefaultpath

clear variables
clc
close all

addpath('./functions')
Powerdbm=20;
Power=db2pow(Powerdbm-30);
f=28*10^9;
c=3*10^8;
lambda=c/f;
k=2*pi/lambda;
Ntile_y=100;
Ntile_z=100;
Qtile_y=1;
Qtile_z=1;
Nant_bs_x=4;
Nant_bs_z=4;
d_bs_x=lambda/2;
d_bs_z=lambda/2;
d_irs_y=lambda/2;
d_irs_z=lambda/2;
Ly=Ntile_y*Qtile_y*lambda/2;
Lz=Ntile_z*Qtile_z*lambda/2;
D=sqrt(Ly^2+Lz^2);
dff=2*D^2/lambda;
dqnf=0.62*sqrt(D^3/lambda);


Rscatter=0;
Lcls=5;
Lsubpath=20;
eta=[2 2 2];
linf=1;
linr=1;
lint=1;


Ground=-6;

p_bs=[30 80 5];
p_irs=[0 0 0];
p_mu=[30 -5 -5];

%pp_bs=[30 80 5;30 80 5];
%pp_mu=[5 5 -5;5 5 -5];
pp_bs=[p_bs;p_bs];
pp_mu=[p_mu;p_mu];

p_mu_virtual=p_mu-[0,0,2*abs(p_mu(1,3)-Ground)];
%errort=180/pi*abs(atan((pt(1,1)-pirsx_center)./abs(pirsy-pt(2,1)))-atan((pt(1,2)-pirsx_center)./abs(pirsy-pt(2,2))))
%errorr=180/pi*abs(atan((pr(1,1)-pirsx_center)./abs(pirsy-pr(2,1)))-atan((pr(1,2)-pirsx_center)./abs(pirsy-pr(2,2))))

%% Ground Reflection
%A=3; B=0; C=0.003; D=0.34; %Very dry ground
A=30.4; B=-0.47; C=0.18; D=1.05; %Medium dry ground
epsilon_r=A*(f/10^9)^B; sigma=C*(f/10^9)^D;
epsilon=epsilon_r-1i*17.98*sigma/(f/10^9);
d_ground=sqrt((p_mu(1,1)-p_irs(1,1))^2+(p_mu(1,2)-p_irs(1,2))^2);
theta_r=atan((p_irs(1,3)-Ground+p_mu(1,3)-Ground)/(d_ground));
Z=sqrt(epsilon-cos(theta_r)^2);
R1=(epsilon*sin(theta_r)-Z)/(epsilon*sin(theta_r)+Z);
R2=(sin(theta_r)-Z)/(sin(theta_r)+Z);
R_ground=sqrt(0.5*abs(R1)^2+0.5*abs(R2)^2);
R_groundfactordb=pow2db(R_ground^2/2);
T_ground=(sqrt((p_irs(1,3)-Ground+p_mu(1,3)-Ground)^2+d_ground^2)-sqrt((p_irs(1,3)-p_mu(1,3))^2+d_ground^2))/c;

%%


Vscatter_d=[0,0,0;40,60,10];
Vscatter_t=[0 0 0;40 50	10];
Vscatter_r=[5 60 5;5 60 5];
Vscatter{1}=Vscatter_d;
Vscatter{2}=Vscatter_t;
Vscatter{3}=Vscatter_r;
for l=1:Lcls
    Pscatter(l,:)=rand(1,3).*[40 10 5]+[0 -10 -10];
end
%Pscatter=[2 51 4;2 49 3];
%% Normalization factors:

h_blk=db2pow(-40);
beta=((lambda/(4*pi))^2);
%beta=db2pow(-46);
d0=1;
sigma2_SF=[0 0 0];
N0=db2pow(-174)/1000;
W=20*10^6;
Nf=db2pow(6);
var_noise_base=N0*W*Nf;
%var_noise_base=3.16978638492224e-13;
norm_factor = var_noise_base;
var_noise = var_noise_base/norm_factor;
factor_Hi = 1/sqrt(sqrt(norm_factor));
factor_Hr = factor_Hi;
factor_Hd = 1/sqrt(norm_factor);

%% Codebook
n=1;
for yy=1:Ntile_y
    for zz=1:Ntile_z
        ppirs(n,:)=[p_irs(1,1) p_irs(1,2)-Ntile_y/2*lambda/2-lambda/4+lambda/2*yy p_irs(1,3)-Ntile_z/2*lambda/2-lambda/4+lambda/2*zz];
        n=n+1;
    end
end

pp_mu_reflector=[Pscatter;p_mu_virtual;pp_mu(1,:)];
for l=1:Lcls+2
temp=[pp_mu_reflector(l,:);pp_mu_reflector(l,:)];
Wna2{l} = near_ana_2D(f,Ntile_y,Ntile_z,temp,pp_bs,p_irs,ppirs);
%Wna2{l}=Wna2{l}*pi; %unit-cell factor=pi
wna2{l}=diag(Wna2{l}); 
end

%% Krice
Krice_dB=[-100,100,0];
sum_reflector(1:5,1:5)=0;
sum_point=0;
sum_fullCSI=0;
sum_specular=0;
sum_random=0;
number_reflector(5,1)=0;



ParamC = struct('p_bs',p_bs,'p_irs',p_irs,'p_mu',p_mu,...
            'factor_Hd',factor_Hd,'factor_Hi',factor_Hi,'factor_Hr',factor_Hr,...
            'eta',eta,'sigma2_SF',sigma2_SF,'d0',d0,'beta',beta,'h_blk',h_blk,...
            'lambda',lambda,'Nant_bs_x',Nant_bs_x,'Nant_bs_z',Nant_bs_z,...
            'Ntile_y',Ntile_y,'Ntile_z',Ntile_z,'Qtile_y',Qtile_y,'Qtile_z',Qtile_z,...
            'd_bs_x',d_bs_x,'d_bs_z',d_bs_z,'d_irs_y',d_irs_y,'d_irs_z',d_irs_z,...
            'Krice_dB',Krice_dB,'Lcls',Lcls,'Lsubpath',Lsubpath,'Rscatter',Rscatter);

        [H_d,H_t,H_r,H_r_los,H_r_nlos,H_r_ground,var_ground] = func_irs_channel_near_modified(ParamC,Vscatter,Pscatter,R_groundfactordb,p_mu_virtual);

hbar_pl_ground=db2pow(R_groundfactordb+var_ground);
ITER=100;
self_blockage=0;
for iter=1:ITER
%% Choose number of scatters and which ones
number_p=floor((ITER-iter)/ITER*5+1);
number_total=1:Lcls;
number_combination=nchoosek(number_total,number_p);
number_random=floor(size(number_combination,1)*rand+1);
number_scatters=number_combination(number_random,:);

Lcls_new=number_p;
Pscatter_new=Pscatter(number_scatters',:);
self_blockage=1-self_blockage;



%% Channel coefficient
for kkk=-8:4:8
   Krice_dB=[-100,100,kkk];
 H_r{1,1} = self_blockage*H_r_los{1,1}+sqrt(1/(db2pow(kkk))-hbar_pl_ground)*H_r_nlos{1,1}+H_r_ground{1,1}; %Power LOS is constant

%% BS precoder and W full CSI
vbs=H_t(1,:)'/norm(H_t(1,:)',2)*sqrt(Power)*pi; %unit-cell
h_t=H_t*vbs;
h_d=H_d{1,1}*vbs;
omega=angle(h_d)*ones(1,Ntile_z*Ntile_y)-angle(h_t).'-angle(H_r{1,1});
w_fullCSI=exp(1i*omega');
W_fullCSI=diag(w_fullCSI);

SNR_fullCSI(kkk/4+3)=(abs((H_d{1,1}+H_r{1,1}*W_fullCSI*H_t)*vbs)^2);
 for pre=1:3
     vbs=(H_d{1,1}+H_r{1,1}*W_fullCSI*H_t)'/norm((H_d{1,1}+H_r{1,1}*W_fullCSI*H_t)',2)*sqrt(Power)*pi; %unit-cell
     omega=angle(h_d)*ones(1,Ntile_z*Ntile_y)-angle(h_t).'-angle(H_r{1,1});
     w_fullCSI=exp(1i*omega');
     W_fullCSI=diag(w_fullCSI);
 end

%% Point focus
 SNR_na2_point(kkk/4+3)=(abs((H_d{1,1}+H_r{1,1}*Wna2{Lcls+2}*H_t)*vbs)^2);
 SNR_thr=SNR_na2_point(kkk/4+3);
 SNR_na2_reflector(kkk/4+3)=SNR_na2_point(kkk/4+3);

%% Reflector focus
SNR_temp=(abs((H_d{1,1}+H_r{1,1}*Wna2{Lcls+1}*H_t)*vbs)^2);
if SNR_temp>SNR_thr
     SNR_na2_reflector(kkk/4+3)=SNR_temp;
     SNR_thr=SNR_temp;
 end
for l=1:Lcls_new
 SNR_temp=(abs((H_d{1,1}+H_r{1,1}*Wna2{number_scatters(l)}*H_t)*vbs)^2);
 if SNR_temp>SNR_thr
     SNR_na2_reflector(kkk/4+3)=SNR_temp;
     SNR_thr=SNR_temp;
 end
end

%% Specular and random
SNR_specular(kkk/4+3)=(abs((H_d{1,1}+H_r{1,1}*diag(ones(Ntile_y*Ntile_z,1))*H_t)*vbs)^2);
SNR_random(kkk/4+3)=(abs((H_d{1,1}+H_r{1,1}*diag(exp(1i*2*pi*rand(Ntile_y*Ntile_z,1)))*H_t)*vbs)^2);


end
number_reflector(number_p)=number_reflector(number_p)+1;
sum_fullCSI=sum_fullCSI+SNR_fullCSI;
sum_point=sum_point+SNR_na2_point;
sum_reflector(number_p,:)=sum_reflector(number_p,:)+SNR_na2_reflector;
sum_specular=sum_specular+SNR_specular;
sum_random=sum_random+SNR_random;
end
sum_fullCSI=pow2db(sum_fullCSI/iter);
sum_point=pow2db(sum_point/iter);
sum_reflector=pow2db(sum_reflector./number_reflector);
sum_specular=pow2db(sum_specular/iter);
sum_random=pow2db(sum_random/iter);
AA=[(-8:4:8)' sum_fullCSI' sum_point' sum_reflector' sum_specular' sum_random'];
kk=(-8:4:8);
plot(kk, sum_fullCSI, kk ,sum_point,kk,sum_reflector(1,:),kk,sum_reflector(2,:),kk,sum_reflector(3,:),kk,sum_reflector(4,:),kk,sum_reflector(5,:),kk,sum_specular,kk,sum_random)


