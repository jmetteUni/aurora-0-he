% This program is used to create initial condition files from NEMO data
%created by G.Xu
%modified by J. Mette, 20/04/2024

clc;
clear;
close all;

%--------------------------------------------------------------------------
%  Initialization
%--------------------------------------------------------------------------
%addpath('../roms_matlab/initial');
%addpath('../roms_matlab/grid');
%addpath('../roms_matlab/utility');
%addpath('../roms_matlab/boundary');
%addpath('../roms_matlab/netcdf');
%addpath('../roms_matlab/mexcdf/mexnc');


app_dir = '/home/jonathan/Dokumente/model/roms_project/aurora-0/';
%GRD_file = 'grid-M512L512.nc';
GRD_file = 'grid-M265L265.nc';
GRD_file = 'grid-M128L128.nc';

GRDname = fullfile(app_dir,'/Data',GRD_file);
NEMO_dir = 'Dokumente/model/inputs/boundary_conditions';
NEMO_name = 'cmems_mod_glo_phy_myint_0.083deg_P1M-m_1719414538489.nc';
NEMOfile = fullfile(NEMO_dir,NEMO_name);


% set output file directory and name
% INI_dir = '/media/jonathan/Extreme SSD/inputs';

% set reference date
time_nemo = ncread(NEMOfile,'time');    %hours since 1950-01-01
Pdate = time_nemo/24+datenum(1950,1,1,0,0,0);
date_ref = Pdate;
date_start = datenum(2022,07,26,00,00,00);

% set output file name
INI_file = sprintf('%s_ini_nemo_%s_64layer.nc',GRD_file(1:end-3),datestr(date_start,'yyyymmddTHH'));
INIname = fullfile(app_dir,'/Data',INI_file);


CREATE = true;                   % logical switch to create NetCDF
report = false;                  % report vertical grid information

% Get number of grid points
Sinp.N = 40;
Sinp.Vtransform = 2;
Sinp.Vstretching = 4;
Sinp.theta_s=0;
Sinp.theta_b=3;
Sinp.Tcline = 20;
Sinp.hc = 20;
G = get_roms_grid(GRDname,Sinp);


%--------------------------------------------------------------------------
%  Set application parameters in structure array, S.
%--------------------------------------------------------------------------
S = G;
S.ncname      = INIname;     % output NetCDF file
S.NT = 2;


%--------------------------------------------------------------------------
%  Set initial conditions for plume tracers
%--------------------------------------------------------------------------

dye_01 = zeros(size(S.z_r));

%--------------------------------------------------------------------------
%  Interpolate initial conditions from NEMO data to application grid.
%--------------------------------------------------------------------------

disp(' ')
disp('Interpolating from NEMO to ROMS grid ...');
disp(' ')


%  Get NEMO grid and fields
Pz = -ncread(NEMOfile,'depth');
Plat = ncread(NEMOfile,'latitude');
Plon = ncread(NEMOfile,'longitude');
Temp = ncread(NEMOfile,'thetao',[1 1 1 1],[inf,inf,inf,1]);
Salt = ncread(NEMOfile,'so',[1 1 1 1],[inf,inf,inf,1]);
Btemp = ncread(NEMOfile,'bottomT',[1 1 1],[inf,inf,1]);
Zeta = ncread(NEMOfile,'zos',[1 1 1],[inf,inf,1]);
Uvel = ncread(NEMOfile,'uo',[1 1 1 1],[inf,inf,inf,1]);
Vvel = ncread(NEMOfile,'vo',[1 1 1 1],[inf,inf,inf,1]);

Plon = double(Plon);
Plat = double(Plat);
Pz = double(Pz);
Zeta = double(Zeta);
Temp = double(Temp);
Btemp = double(Btemp);
Salt = double(Salt);
Uvel = double(Uvel);
Vvel = double(Vvel);
Zeta_p = permute(Zeta,[2,1]);
Temp_p = permute(Temp,[2,1,3]);
Btemp_p = permute(Btemp,[2,1]);
Salt_p = permute(Salt,[2,1,3]);
Uvel_p = permute(Uvel,[2,1,3]);
Vvel_p = permute(Vvel,[2,1,3]);
clear Zeta Temp Btemp Salt Uvel Vvel

%  Set initial conditions time (seconds).
S.ocean_time = (Pdate-date_ref)*24*3600;

%  Interpolate free-surface initial conditions.
Tlon = S.lon_rho;
Tlat = S.lat_rho;
[Pllon,Pllat] = meshgrid(Plon,Plat);
dummy = ones(size(Pllon));
F_2D = scatteredInterpolant(Pllon(:),Pllat(:),dummy(:));
F_2D.Values = Zeta_p(:);
zeta = F_2D(Tlon,Tlat);
F_2D.Values = Btemp_p(:);
btemp = F_2D(Tlon,Tlat);

%  Interpolate 3D variables
zind = find(Pz<min(S.z_r(:)),1);
Temp_int = zeros(size(Tlon,1),size(Tlon,2),zind);
Salt_int = zeros(size(Tlon,1),size(Tlon,2),zind);
U_int = zeros(size(Tlon,1),size(Tlon,2),zind);
V_int = zeros(size(Tlon,1),size(Tlon,2),zind);
for e = 1:zind

    % first round of interpolation
    Temp1 = Temp_p(:,:,e);
    Salt1 = Salt_p(:,:,e);
    Uvel1 = Uvel_p(:,:,e);
    Vvel1 = Vvel_p(:,:,e);
    F_2D.Values = Temp1(:);
    Temp_int1 = F_2D(Tlon,Tlat);
    F_2D.Values = Salt1(:);
    Salt_int1 = F_2D(Tlon,Tlat);
    F_2D.Values = Uvel1(:);
    U_int1 = F_2D(Tlon,Tlat);
    F_2D.Values = Vvel1(:);
    V_int1 = F_2D(Tlon,Tlat);

    %second round to replace NaNs with the nearest valid values
        % ind = find(isnan(Temp_int1(:)));
        % if any(isnan(Temp_int1(:))) && ~all(isnan(Temp_int1(:)))
        %    ii = find(~isnan(Temp1(:)));
        %    F_2D2 = scatteredInterpolant(Pllon(ii),Pllat(ii),Temp1(ii),'nearest');
        %    Temp_int1(ind) = F_2D2(Tlon(ind),Tlat(ind));
        %    F_2D2.Values = Salt1(ii);
        %    Salt_int1(ind) = F_2D2(Tlon(ind),Tlat(ind));
        %    ii = find(~isnan(Uvel1(:)));
        %    F_2D2 = scatteredInterpolant(Pllon(ii),Pllat(ii),Uvel1(ii),'nearest');
        %    U_int1(ind) = F_2D2(Tlon(ind),Tlat(ind));
        %    F_2D2.Values = Vvel1(ii);
        %    V_int1(ind) = F_2D2(Tlon(ind),Tlat(ind));
        % end

    Temp_int(:,:,e) = Temp_int1;
    Salt_int(:,:,e) = Salt_int1;
    U_int(:,:,e) = U_int1;
    V_int(:,:,e) = V_int1;
end


% for e = 1:zind
%     Temp2 = Temp_int1(:,:,e);
%     ind = find(isnan(Temp2));
%     if ~isempty(ind)
%        ii = find(~isnan(Temp1(:)));
%        F_2D2 = scatteredInterpolant(Pllon(ii),Pllat(ii),Temp1(ii),'nearest');
%        Temp_int1(:,:,e) = 
%     end
%     ind2 = find(~isnan(Temp2));
%     F_2D2 = scatteredInterpolant(Tlon(ind2),Tlat(ind2),Temp2(ind2),'nearest');
%     Temp2(ind) = F_2D2(Tlon(ind),Tlat(ind));
%     Temp_int1(:,:,e) = Temp2;
%
%     Salt2 = Salt_int1(:,:,e);
%     F_2D2.Values = Salt2(ind2);
%     Salt2(ind) = F_2D2(Tlon(ind),Tlat(ind));
%     Salt_int1(:,:,e) = Salt2;
%
%     Uvel2 = U_int1(:,:,e);
%     F_2D2.Values = Uvel2(ind2);
%     Uvel2(ind) = F_2D2(Tlon(ind),Tlat(ind));
%     U_int1(:,:,e) = Uvel2;
%
%     Vvel2 = V_int1(:,:,e);
%     F_2D2.Values = Vvel2(ind2);
%     Vvel2(ind) = F_2D2(Tlon(ind),Tlat(ind));
%     V_int1(:,:,e) = Vvel2;
% end

temp = zeros(size(S.z_r));
salt = zeros(size(S.z_r));
Urho = zeros(size(S.z_r));
Vrho = zeros(size(S.z_r));
for i = 1:size(Tlon,1)
    for j = 1:size(Tlon,2)
        h1 = -G.h(i,j);
        z_r1 = squeeze(S.z_r(i,j,:));
        temp1 = squeeze(Temp_int(i,j,:));
        salt1 = squeeze(Salt_int(i,j,:));
        U1 = squeeze(U_int(i,j,:));
        V1 = squeeze(V_int(i,j,:));
        ikeep_ts = find(~isnan(temp1));
        ikeep_uv = find(~isnan(U1));
        depth_int1_ts = [Pz(ikeep_ts);h1];
        depth_int1_uv = [Pz(ikeep_uv);h1];
        Temp_int1 = [temp1(ikeep_ts);btemp(i,j)];
        Salt_int1 = [salt1(ikeep_ts);salt1(ikeep_ts(end))];
        U_int1 = [U1(ikeep_uv);0];
        V_int1 = [V1(ikeep_uv);0];
        %         temp(i,j,:) = interp1(Pdepth(1:zind),squeeze(Temp_int(i,j,:)),squeeze(S.z_r(i,j,:)));
        %         salt(i,j,:) = interp1(Pdepth(1:zind),squeeze(Salt_int(i,j,:)),squeeze(S.z_r(i,j,:)));
        %         Urho(i,j,:) = interp1(Pdepth(1:zind),squeeze(U_int(i,j,:)),squeeze(S.z_r(i,j,:)));
        %         Vrho(i,j,:) = interp1(Pdepth(1:zind),squeeze(V_int(i,j,:)),squeeze(S.z_r(i,j,:)));
        temp(i,j,:) = interp1(depth_int1_ts(:),Temp_int1(:),z_r1(:));
        salt(i,j,:) = interp1(depth_int1_ts(:),Salt_int1(:),z_r1(:));
        Urho(i,j,:) = interp1(depth_int1_uv(:),U_int1(:),z_r1(:));
        Vrho(i,j,:) = interp1(depth_int1_uv(:),V_int1(:),z_r1(:));
    end
end


%  Process velocity: rotate and/or average to staggered C-grid locations.

[u,v]=roms_vectors(Urho,Vrho,S.angle,S.mask_u,S.mask_v);

%  Compute barotropic velocities by vertically integrating (u,v).

[ubar,vbar]=uv_barotropic(u,v,S.Hz);

% check if any variables have NaNs, if yes interpolate
inan_temp = find(isnan(temp(:)));
if ~isempty(inan_temp)
  nantemp = isnan(temp);
  t = 1:numel(temp);
  temp(nantemp) = interp1(t(~nantemp), temp(~nantemp), t(nantemp));
  fprintf('NaNs found and interpolated in Temp');
end
inan_salt = find(isnan(salt(:)));
if ~isempty(inan_salt)
  nansalt = isnan(salt);
  t = 1:numel(salt);
  salt(nansalt) = interp1(t(~nansalt), salt(~nansalt), t(nansalt));
  fprintf('NaNs found and interpolated in Salt');
end
inan_u = find(isnan(u(:)));
if ~isempty(inan_u)
  nanu = isnan(u);
  t = 1:numel(u);
  u(nanu) = interp1(t(~nanu), u(~nanu), t(nanu));
  fprintf('NaNs found and interpolated in u');
end
inan_v = find(isnan(v(:)));
if ~isempty(inan_v)
  nanv = isnan(v);
  t = 1:numel(v);
  v(nanv) = interp1(t(~nanv), v(~nanv), t(nanv));
  fprintf('NaNs found and interpolated in v');
end
inan_zeta = find(isnan(zeta(:)));
if ~isempty(inan_zeta)
  nanzeta = isnan(zeta);
  t = 1:numel(zeta);
  zeta(nanzeta) = interp1(t(~nanzeta), zeta(~nanzeta), t(nanzeta));
  fprintf('NaNs found and interpolated in Zeta');
end

%--------------------------------------------------------------------------
%  Create initial condition Netcdf file.
%--------------------------------------------------------------------------

if (CREATE)
    [status]=c_initial(S);

    %  Set attributes for "ocean_time".

    avalue='seconds since 2022-07-26 00:00:00';
    %avalue = sprintf('seconds since %s',datestr(date_ref,'yyyy-mm-dd HH:MM:SS'));
    [status]=nc_attadd(INIname,'units',avalue,'ocean_time');

    avalue='gregorian';
    [status]=nc_attadd(INIname,'calendar',avalue,'ocean_time');

    %  Set global attributes.

    avalue='Axial Seamount';
    [status]=nc_attadd(INIname,'title',avalue);

    avalue='global-reanalysis-phy-001-030-daily';
    [status]=nc_attadd(INIname,'data_source',avalue);

    [status]=nc_attadd(INIname,'grd_file',GRDname);
end

%--------------------------------------------------------------------------
%  Write out initial conditions.
%--------------------------------------------------------------------------

if (CREATE)
    disp(' ')
    disp([ 'Writing initial conditions ...']);
    disp(' ')

    [status]=nc_write(INIname, 'spherical',   S.spherical);
    [status]=nc_write(INIname, 'Vtransform',  S.Vtransform);
    [status]=nc_write(INIname, 'Vstretching', S.Vstretching);
    [status]=nc_write(INIname, 'theta_s',     S.theta_s);
    [status]=nc_write(INIname, 'theta_b',     S.theta_b);
    [status]=nc_write(INIname, 'Tcline',      S.Tcline);
    [status]=nc_write(INIname, 'hc',          S.hc);
    [status]=nc_write(INIname, 's_rho',       S.s_rho);
    [status]=nc_write(INIname, 's_w',         S.s_w);
    [status]=nc_write(INIname, 'Cs_r',        S.Cs_r);
    [status]=nc_write(INIname, 'Cs_w',        S.Cs_w);

    [status]=nc_write(INIname, 'h',           S.h);
    [status]=nc_write(INIname, 'lon_rho',     S.lon_rho);
    [status]=nc_write(INIname, 'lat_rho',     S.lat_rho);
    [status]=nc_write(INIname, 'lon_u',       S.lon_u);
    [status]=nc_write(INIname, 'lat_u',       S.lat_u);
    [status]=nc_write(INIname, 'lon_v',       S.lon_v);
    [status]=nc_write(INIname, 'lat_v',       S.lat_v);

    IniRec = 1;

    [status]=nc_write(INIname, 'ocean_time', S.ocean_time, IniRec);

    [status]=nc_write(INIname, 'zeta', zeta, IniRec);
    [status]=nc_write(INIname, 'ubar', ubar, IniRec);
    [status]=nc_write(INIname, 'vbar', vbar, IniRec);
    [status]=nc_write(INIname, 'u',    u,    IniRec);
    [status]=nc_write(INIname, 'v',    v,    IniRec);
    [status]=nc_write(INIname, 'temp', temp, IniRec);
    [status]=nc_write(INIname, 'salt', salt, IniRec);
    [status]=nc_write(INIname, 'dye_01', dye_01, IniRec);

end


