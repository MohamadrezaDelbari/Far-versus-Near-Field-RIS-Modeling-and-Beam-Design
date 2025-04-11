%%----------------------------------------------------%%
%%----- Mohamadreza Delbari
%%      Please cite our paper:
%%----- DOI: https://arxiv.org/pdf/2401.08237
%%----------------------------------------------------%%
function [pp_bs_ant,pp_irs_uc_vec,pp_tile,pp_irs_uc_tile] = func_bs_irs_ant_p(Param)

Nant_bs_x = Param.Nant_bs_x; % number of BS's antenna along x axis
Nant_bs_z = Param.Nant_bs_z; % number of BS's antenna along z axis
Nant_bs = Nant_bs_x*Nant_bs_z;

d_bs_x = Param.d_bs_x; % BS element spacing along x axis
d_bs_z = Param.d_bs_z; % BS element spacing along z axis

Ntile_y = Param.Ntile_y; % number of IRS tiles along y axis
Ntile_z = Param.Ntile_z; % number of IRS tiles along y axis
Ntile = Ntile_y*Ntile_z;

Qtile_y = Param.Qtile_y; % number of elements per tile along y axis
Qtile_z = Param.Qtile_z; % number of elements per tile along y axis
Qtile = Qtile_y*Qtile_z;

d_irs_y = Param.d_irs_y; % IRS element spacing along y axis
d_irs_z = Param.d_irs_z; % IRS element spacing along y axis

lambda = Param.lambda; % wavelength

%% Antenna positions:
% positions are computed in the respective coordinate system of the nodes

% BS antennas
pp_bs_ant_x = linspace(-(Nant_bs_x-1)*d_bs_x/2,(Nant_bs_x-1)*d_bs_x/2,Nant_bs_x);
pp_bs_ant_z = linspace(-(Nant_bs_z-1)*d_bs_z/2,(Nant_bs_z-1)*d_bs_z/2,Nant_bs_z);
[XX,ZZ] = meshgrid(pp_bs_ant_x,pp_bs_ant_z);
vx=reshape(XX,[],1);
vz=reshape(ZZ,[],1);
pp_bs_ant = [vx,zeros(Nant_bs,1),vz];

% for nn=1:Nant_bs
%         plot(pp_bs_ant(nn,1)/d_bs_x,pp_bs_ant(nn,3)/d_bs_z,'o')
%         hold on
% end

% IRS tiles
d_tile_y = Qtile_y*d_irs_y; % tile spacing
d_tile_z = Qtile_z*d_irs_z; % tile spacing

pp_tile_y = linspace(-(Ntile_y-1)*d_tile_y/2,(Ntile_y-1)*d_tile_y/2,Ntile_y);
pp_tile_z = linspace(-(Ntile_z-1)*d_tile_z/2,(Ntile_z-1)*d_tile_z/2,Ntile_z);
[YY,ZZ] = meshgrid(pp_tile_y,pp_tile_z);
vy=reshape(YY,[],1);
vz=reshape(ZZ,[],1);
pp_tile = [zeros(Ntile,1),vy,vz];

% for nn=1:Ntile
%         plot(pp_tile(nn,2)/d_tile_y,pp_tile(nn,3)/d_tile_z,'o')
%         hold on
% end

% Tile unit-cells
pp_ref_y = linspace(-(Qtile_y-1)*d_irs_y/2,(Qtile_y-1)*d_irs_y/2,Qtile_y);
pp_ref_z = linspace(-(Qtile_z-1)*d_irs_z/2,(Qtile_z-1)*d_irs_z/2,Qtile_z);
[YY,ZZ] = meshgrid(pp_ref_y,pp_ref_z);
vy=reshape(YY,[],1);
vz=reshape(ZZ,[],1);
pp_ref = [zeros(Qtile,1),vy,vz];

% for nn=1:Qtile
%         plot(pp_ref(nn,2)/d_irs_y,pp_ref(nn,3)/d_irs_z,'o')
%         hold on
% end

pp_irs_uc_vec = [];
for tt=1:Ntile
    pp_irs_uc_tile{tt} = pp_tile(tt,:) + pp_ref;
    pp_irs_uc_vec = [pp_irs_uc_vec; pp_irs_uc_tile{tt}];
end

% for tt=1:Ntile
%     plot(pp_irs_uc{tt}(:,2)/d_irs_y,pp_irs_uc{tt}(:,3)/d_irs_z,'o')
%         hold on
%     plot(pp_tile(tt,2)/d_irs_y,pp_tile(tt,3)/d_irs_z,'x')
%         hold on
% end

end
