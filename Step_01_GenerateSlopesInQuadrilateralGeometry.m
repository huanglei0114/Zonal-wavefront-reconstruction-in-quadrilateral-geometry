% Generate Slopes in Quadrilateral Geometry.

%% Clear everything.

close all
clear
clc

%% Define analytical height and slopes.

syms x y

z = 0.1*...
   ( 3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) ... 
   - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) ... 
   - 1/3*exp(-(x+1).^2 - y.^2) ...
   );
sx = diff(z,x);
sy = diff(z,y);

% size
Nx = 128;
Ny = 128;

% range in x and y
Min = -2;
Max = +2;
[x, y] = meshgrid(linspace(Min,Max,Nx),linspace(Min,Max,Ny));
x0 = x;
y0 = y;

% random position errors.
dy = diff(y0);
dx = diff(x0,[],2);
std_pos = 5e-2*median([dy(:);dx(:)]);
rpex = std_pos*randn(size(x0));
rpey = std_pos*randn(size(y0));

% distortion
[~, r] = cart2pol(x0,y0);
kd = -1e-2;
dx = kd.*r.^2.*x0;
dy = kd.*r.^2.*y0;

% keystone
theta_x = -20;
yt = -4;

R = [1              0               0
     0              cosd(theta_x)   -sind(theta_x)
     0              sind(theta_x)   cosd(theta_x)];
 
X0 = [x0(:), y0(:), x0(:)*0]';
Xr = R*X0;
xr = reshape(Xr(1,:),Ny,Nx);
yr = reshape(Xr(2,:),Ny,Nx);
zr = reshape(Xr(3,:),Ny,Nx);

T = [0; yt; -yt./tand(theta_x)];
xk = (xr-T(1))./(zr-T(3))*(-T(3))+T(1);
yk = (yr-T(2))./(zr-T(3))*(-T(3))+T(2);

kx = xk - x0;
ky = yk - y0;

% implement 
x = x0 + rpex + dx + kx;
y = y0 + rpey + dy + ky;

% evaluate the height and slopes
z_anal  = eval( z);
sx_anal = eval(sx);
sy_anal = eval(sy);

% Invalid Regions, if want to simulate incomplete dataset.
% NanMask = sx_anal'>0.5;
% z_anal(NanMask) = NaN;
% sx_anal(NanMask) = NaN;
% sy_anal(NanMask) = NaN;

FontSize = 20;
figure('Name','Nominal z'); 
surf(x,y,z_anal);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
title('z [mm]');
view([0,90]);
axis image;
box on;
shading interp;
set(gca,'FontSize',FontSize,'FontWeight','Bold')
colorbar;

figure('Name','Nominal x-slope'); 
surf(x,y,sx_anal);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('x-slope [rad]');
title('x-slope [rad]');
view([0,90]);
axis image;
box on;
shading interp;
set(gca,'FontSize',FontSize,'FontWeight','Bold')
colorbar;
colormap(RGBColorMap);

figure('Name','Nominal y-slope'); 
surf(x,y,sy_anal);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('y-slope [rad]');
title('y-slope [rad]');
view([0,90]);
axis image;
box on;
shading interp;
set(gca,'FontSize',FontSize,'FontWeight','Bold')
colorbar;
colormap(RGBColorMap);


save Step_01_GenerateSlopesInQuadrilateralGeometry.mat ...
    x y sx_anal sy_anal z_anal
