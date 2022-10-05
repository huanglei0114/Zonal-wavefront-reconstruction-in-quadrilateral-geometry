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
Nx = 16;
Ny = 16;

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

% implement................................................................ 
x = x0 + rpex;
y = y0 + rpey;

figure; 
mesh(x0,y0,x0*0,'linewidth',2,'facecolor','none','edgecolor',[0.2 0.8 0.2]);
hold on;
mesh(x,y,x*0,'linewidth',1,'facecolor','none','edgecolor',[0.8 0.2 0.8] ...
    , 'marker','.','markerfacecolor',[0.8 0.2 0.8],'markersize',16);
view([0, 90]);
axis image off; 


% implement................................................................ 
x = x0 + dx;
y = y0 + dy;

figure; 
mesh(x0,y0,x0*0,'linewidth',2,'facecolor','none','edgecolor',[0.2 0.8 0.2]);
hold on;
mesh(x,y,x*0,'linewidth',1,'facecolor','none','edgecolor',[0.8 0.2 0.8] ...
    , 'marker','.','markerfacecolor',[0.8 0.2 0.8],'markersize',16);
view([0, 90]);
axis image off; 


% implement................................................................ 
x = x0 + kx;
y = y0 + ky;

figure; 
mesh(x0,y0,x0*0,'linewidth',2,'facecolor','none','edgecolor',[0.2 0.8 0.2]);
hold on;
mesh(x,y,x*0,'linewidth',1,'facecolor','none','edgecolor',[0.8 0.2 0.8] ...
    , 'marker','.','markerfacecolor',[0.8 0.2 0.8],'markersize',16);
view([0, 90]);
axis image off; 


% implement................................................................ 
x = x0 + rpex + dx + kx;
y = y0 + rpey + dy + ky;

figure; 
mesh(x0,y0,x0*0,'linewidth',2,'facecolor','none','edgecolor',[0.2 0.8 0.2]);
hold on;
mesh(x,y,x*0,'linewidth',1,'facecolor','none','edgecolor',[0.8 0.2 0.8] ...
    , 'marker','.','markerfacecolor',[0.8 0.2 0.8],'markersize',16);
view([0, 90]);
axis image off; 
