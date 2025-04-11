%%----------------------------------------------------%%
%%----- Mohamadreza Delbari
%%      Please cite our paper:
%%----- DOI: https://arxiv.org/pdf/2401.08237
%%----------------------------------------------------%%
function [H_d,H_i,H_r,H_r_los,H_r_nlos,H_r_ground,var_ground] = func_irs_channel_near_modified(Param,Vscatter,Pscatter,R_groundfactordb,p_mu_virtual)

p_bs = Param.p_bs; % BS's positions
p_irs = Param.p_irs; % IRS's positions
p_mu = Param.p_mu; % MU's positions

Nant_bs_x = Param.Nant_bs_x; % number of BS's antenna along x axis
Nant_bs_z = Param.Nant_bs_z; % number of BS's antenna along z axis
Nant_bs = Nant_bs_x*Nant_bs_z;

d_bs_x = Param.d_bs_x; % BS element spacing along x axis
d_bs_z = Param.d_bs_z; % BS element spacing along z axis

Nant_mu = 1; % number of MU's antenna

Ntile_y = Param.Ntile_y; % number of IRS tiles along y axis
Ntile_z = Param.Ntile_z; % number of IRS tiles along y axis
Ntile = Ntile_y*Ntile_z;

Qtile_y = Param.Qtile_y; % number of elements per tile along y axis
Qtile_z = Param.Qtile_z; % number of elements per tile along y axis
Qtile = Qtile_y*Qtile_z;

d_irs_y = Param.d_irs_y; % IRS element spacing along y axis
d_irs_z = Param.d_irs_z; % IRS element spacing along y axis

Krice_dB = Param.Krice_dB; % K factor 

factor_Hd = Param.factor_Hd; % normalization factor
factor_Hi = Param.factor_Hi; % normalization factor
factor_Hr = Param.factor_Hr; % normalization factor

eta = Param.eta; % pathloss exponent
sigma2_SF = Param.sigma2_SF; % shadow fading variance
d0 = Param.d0; % reference distance
beta = Param.beta; % reference loss
h_blk = Param.h_blk; % blockage gain for dircet link

lambda = Param.lambda; % wavelength

Lcls = Param.Lcls; % number of non-LOS clusters
Lsubpath = Param.Lsubpath; % number of non-LOS paths
Rscatter = Param.Rscatter; % scattering radius for subpaths

[K,~] = size(p_mu); % number of users

eta_d = eta(1); % pathloss exponent for BS-MU link
eta_i = eta(2); % pathloss exponent for BS-IRS link
eta_r = eta(3); % pathloss exponent for IRS-MU link

sigma2_SF_d = sigma2_SF(1); % show fading variance for BS-MU link
sigma2_SF_i = sigma2_SF(2); % show fading variance for BS-IRS  link
sigma2_SF_r = sigma2_SF(3); % show fading variance for IRS-MU link

kappa = 2*pi/lambda; % wavenumber

Krice = db2pow(Krice_dB); % Rician factor
Krice_d = Krice(1); % Rician factor for BS-MU link
Krice_i = Krice(2); % Rician factorfor BS-IRS  link
Krice_r = Krice(3); % Rician factor for IRS-MU link

Vscatter_d = Vscatter{1}; % scattering volume for paths in BS-MU link
Vscatter_i = Vscatter{2}; % scattering volume for paths in BS-MU link
Vscatter_r = Vscatter{3}; % scattering volume for paths in BS-MU link

Lpath_tot = Lcls*Lsubpath; % total number of paths

%% large scale pathloss and fading:
for kk=1:K
    d = norm(p_bs-p_mu(kk,:));
    hbar_pl_d(kk) = beta*(d0/d)^eta_d; % pathloss
    hbar_sf_d(kk) = 10^(0.1*(normrnd(0,sqrt(sigma2_SF_d)))); % shadow fading

    d = norm(p_irs-p_mu(kk,:));
    hbar_pl_r(kk) = beta*(d0/d)^eta_r; % pathloss
    hbar_sf_r(kk) = 10^(0.1*(normrnd(0,sqrt(sigma2_SF_r)))); % shadow fading

    var_ground=randn;
    hbar_pl_ground(kk)=db2pow(R_groundfactordb+var_ground)*hbar_pl_r(kk);
end

d = norm(p_bs-p_irs);
hbar_pl_i = beta*(d0/d)^eta_i; % pathloss
hbar_sf_i = 10^(0.1*(normrnd(0,sqrt(sigma2_SF_i)))); % shadow fading



%% Positions

% antenna positions:
% positions are computed in the respective coordinate system of the nodes
ParamP = struct('lambda',lambda,'Nant_bs_x',Nant_bs_x,'Nant_bs_z',Nant_bs_z,...
        'Ntile_y',Ntile_y,'Ntile_z',Ntile_z,'Qtile_y',Qtile_y,'Qtile_z',Qtile_z,...
        'd_bs_x',d_bs_x,'d_bs_z',d_bs_z,'d_irs_y',d_irs_y,'d_irs_z',d_irs_z);
[pp_bs_ant,pp_irs_uc_vec,pp_tile,pp_irs_uc_tile]  = func_bs_irs_ant_p(ParamP);

pp_tile_abs = pp_tile+p_irs; % position of tile in the global coordinate system
pp_irs_uc_tile_base = pp_irs_uc_tile{1}-pp_tile(1,:); % unit-cell position in the tile local coordinate system

% scattering positions:
for kk=1:K
pp_path_d{kk} = [];
pp_path_i = [];
pp_path_r{kk} = [];

pp_scr_d{kk} = [];
pp_scr_i = [];
pp_scr_r{kk} = [];
end

for ll=1:Lcls
    pp_scr_i = (Vscatter_i(2,:)-Vscatter_i(1,:)).*rand(1,3)+Vscatter_i(1,:); % scattering cluster position
    for kk=1:K
    pp_scr_d{kk} = (Vscatter_d(2,:)-Vscatter_d(1,:)).*rand(1,3)+Vscatter_d(1,:); % scattering cluster position
    %pp_scr_r{kk} = (Vscatter_r(2,:)-Vscatter_r(1,:)).*rand(1,3)+Vscatter_r(1,:); % scattering cluster position
    pp_scr_r{kk} = Pscatter(ll,:);
    end
    for ss=1:Lsubpath
        temp = (rand(1,3)-0.5)*2;
        pp_path_i = [pp_path_i; pp_scr_i+temp*Rscatter]; % scattering subpaths
        for kk=1:K
        temp = (rand(1,3)-0.5)*2;
        pp_path_d{kk} = [pp_path_d{kk}; pp_scr_d{kk}+temp*Rscatter]; % scattering subpaths
        temp = (rand(1,3)-0.5)*2;
        pp_path_r{kk} = [pp_path_r{kk}; pp_scr_r{kk}+temp*Rscatter]; % scattering subpaths
        end
    end    
end




%% Steering vectors:

% for ll=1:Lpath_tot
%    plot3(pp_path_d(ll,1),pp_path_d(ll,2),(pp_path_d(ll,3)),'o')
%    hold on
% end

% BS-IRS steering vector:
a_bs_irs = zeros(Nant_bs,Lpath_tot+1);
d_BStoIRS = (p_irs-p_bs)/norm(p_irs-p_bs);
a_bs_irs(:,1) = exp(1j*kappa*pp_bs_ant*d_BStoIRS');
for ll=1:Lpath_tot
    d_BStoScr = (pp_path_i(ll,:)-p_bs)/norm(pp_path_i(ll,:)-p_bs);
    a_bs_irs(:,ll+1) = exp(1j*kappa*pp_bs_ant*d_BStoScr');
end

% IRS-BS steering vectors:
a_irs_bs_vec = [];
for tt=1:Ntile
    d_TiletoBS = (p_bs-pp_tile_abs(tt))/norm(p_bs-pp_tile_abs(tt)); %% near field
    phase_0 = mod(kappa*norm(p_bs-pp_tile_abs(tt,:)),2*pi);
    a_irs_bs_tile{tt} = exp(1j*kappa*pp_irs_uc_tile_base*d_TiletoBS'+1j*phase_0);
    a_irs_bs_vec = [a_irs_bs_vec;a_irs_bs_tile{tt}];
end
for ll=1:Lpath_tot
    a_irs_bs_vec_temp = [];
    for tt=1:Ntile
        d_TiletoScr = (pp_path_i(ll,:)-pp_tile_abs(tt,:))/norm(pp_path_i(ll,:)-pp_tile_abs(tt));%% near field
        phase_0 = mod(kappa*norm(pp_path_i(ll,:)-pp_tile_abs(tt)),2*pi);
        a_irs_bs_tile{tt} = [a_irs_bs_tile{tt} ...
            exp(1j*kappa*pp_irs_uc_tile_base*d_TiletoScr'+1j*phase_0)];
        a_irs_bs_vec_temp = [a_irs_bs_vec_temp;a_irs_bs_tile{tt}(:,ll+1)];
    end
    a_irs_bs_vec = [a_irs_bs_vec a_irs_bs_vec_temp];
end


% BS-MU steering vectors:
for kk=1:K
    d_BStoMU = (p_mu(kk,:)-p_bs)/norm(p_mu(kk,:)-p_bs);
    a_bs_mu{kk} = exp(1j*kappa*pp_bs_ant*d_BStoMU');
    for ll=1:Lpath_tot
         d_BStoScr = (pp_path_d{kk}(ll,:)-p_bs)/norm(pp_path_d{kk}(ll,:)-p_bs);
         a_bs_mu{kk}(:,ll+1) = exp(1j*kappa*pp_bs_ant*d_BStoScr');
    end
end

% IRS-MU steering vectors:
for kk=1:K
    a_irs_mu_vec{kk} = [];
    for tt=1:Ntile
        d_TiletoMU = (p_mu(kk,:)-pp_tile_abs(tt,:))/norm(p_mu(kk,:)-pp_tile_abs(tt)); %% near field
        phase_0 = mod(kappa*norm(p_mu(kk,:)-pp_tile_abs(tt,:)),2*pi);
        a_irs_mu_tile{tt,kk} = exp(1j*kappa*pp_irs_uc_tile_base*d_TiletoMU'+1j*phase_0);
        a_irs_mu_vec{kk} = [a_irs_mu_vec{kk};a_irs_mu_tile{tt,kk}];
    end

    for ll=1:Lpath_tot
        a_irs_mu_vec_temp = [];
        for tt=1:Ntile
            d_TiletoScr = (pp_path_r{kk}(ll,:)-pp_tile_abs(tt,:))/norm(pp_path_r{kk}(ll,:)-pp_tile_abs(tt)); %% near field
            phase_0 = mod(kappa*norm(pp_path_r{kk}(ll,:)-pp_tile_abs(tt,:)),2*pi);
            a_irs_mu_tile{tt,kk} = [a_irs_mu_tile{tt,kk}... 
                exp(1j*kappa*pp_irs_uc_tile_base*d_TiletoScr'+1j*phase_0)];
            a_irs_mu_vec_temp = [a_irs_mu_vec_temp;a_irs_mu_tile{tt,kk}(:,ll+1)];
        end
        a_irs_mu_vec{kk} = [a_irs_mu_vec{kk} a_irs_mu_vec_temp];
    end
    
    a_irs_ground_temp=[];
    for tt=1:Ntile
        d_Tiletoground=(p_mu_virtual-pp_tile_abs(tt,:))/norm(p_mu_virtual-pp_tile_abs(tt)); %% near field
        phase_0 = mod(kappa*norm(p_mu_virtual-pp_tile_abs(tt,:)),2*pi);
        a_irs_ground_tile(tt)=exp(1i*phase_0);
        a_irs_ground_temp=[a_irs_ground_temp;a_irs_ground_tile(tt)];
    end
    a_irs_ground_vec{kk}=a_irs_ground_temp;
end


% Hi LOS channel:
H_i_los = factor_Hi*sqrt(hbar_pl_i*hbar_sf_i)*a_irs_bs_vec(:,1)*a_bs_irs(:,1)';
Hpath_i = 1/sqrt(Lpath_tot)*diag(exp(1j*2*pi*rand(1,Lpath_tot))); % random phase for paths
H_i_nlos = factor_Hi*sqrt(hbar_pl_i*hbar_sf_i)*a_irs_bs_vec(:,2:end)*Hpath_i*a_bs_irs(:,2:end)';
for tt=1:Ntile
    H_i_los_tile{tt} = factor_Hi*sqrt(hbar_pl_i*hbar_sf_i)*a_irs_bs_tile{tt}(:,1)*a_bs_irs(:,1)';
    H_i_nlos_tile{tt} = factor_Hi*sqrt(hbar_pl_i*hbar_sf_i)*a_irs_bs_tile{tt}(:,2:end)*Hpath_i*a_bs_irs(:,2:end)';
end


% Hr LOS channel:
for kk=1:K
    H_r_los{kk} = factor_Hr*sqrt(hbar_pl_r(kk)*hbar_sf_r(kk))*a_irs_mu_vec{kk}(:,1)';
    Hpath_r = 1/sqrt(Lpath_tot)*exp(1j*2*pi*rand(1,Lpath_tot)); % random phase for paths
    %Hpath_r = exp(1j*2*pi*rand(1,Lpath_tot));
    H_r_nlos{kk} = factor_Hr*sqrt(hbar_pl_r(kk)*hbar_sf_r(kk))*Hpath_r*a_irs_mu_vec{kk}(:,2:end)';
    H_r_ground{kk}=factor_Hr*sqrt(hbar_pl_ground(kk))*a_irs_ground_vec{kk}(:,1)';
    for tt=1:Ntile
        H_r_los_tile{tt,kk} = factor_Hr*sqrt(hbar_pl_r(kk)*hbar_sf_r(kk))*a_irs_mu_tile{tt,kk}(:,1)';
        H_r_nlos_tile{tt,kk} = factor_Hr*sqrt(hbar_pl_r(kk)*hbar_sf_r(kk))*Hpath_r*a_irs_mu_tile{tt,kk}(:,2:end)';
    end
end

% Hd LOS channel:
for kk=1:K
    H_d_los{kk} = factor_Hd*sqrt(h_blk*hbar_pl_d(kk)*hbar_sf_d(kk))*a_bs_mu{kk}(:,1)';
    Hpath_d = 1/sqrt(Lpath_tot)*exp(1j*2*pi*rand(1,Lpath_tot)); % random phase for paths
    H_d_nlos{kk} = factor_Hd*sqrt(h_blk*hbar_pl_d(kk)*hbar_sf_d(kk))*Hpath_d*a_bs_mu{kk}(:,2:end)';
end

%% Combined channel:
H_i = sqrt(Krice_i/(1+Krice_i))*H_i_los+sqrt(1/(1+Krice_i))*H_i_nlos;

for kk=1:K
    H_d{kk} = sqrt(Krice_d/(1+Krice_d))*H_d_los{kk}+sqrt(1/(1+Krice_d))*H_d_nlos{kk};
    %H_r{kk} =
    %sqrt(Krice_r/(1+Krice_r))*H_r_los{kk}+sqrt(1/(1+Krice_r))*H_r_nlos{kk};
    %total power is constant
    H_r{kk} = floor(rand*2)*H_r_los{kk}+sqrt(1/(Krice_r))*H_r_nlos{kk}+H_r_ground{kk}; %Power LOS is constant
end
