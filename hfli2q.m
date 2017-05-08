function z_hfli2q = hfli2q(sx, sy, x, y, z)
%HFLI2Q Higher-order Finite-difference-based Least-squares Integration in 
%   quadrilateral geometry. (Algorithm 1 in Reference)
%   D * Z = G.
%
%   Reference: 
%   Guanghui Li, Yanqiu Li, Ke Liu, Xu Ma, and Hai Wang, "Improving 
%   wavefront reconstruction accuracy by using integration equations with 
%   higher-order truncation errors in the Southwell geometry," J. Opt. Soc.
%   Am. A 30, 1448-1459 (2013) 

%   Copyright since 2016 by Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%   2016-10-01 Original Version

% Check the number of arguments............................................
% Validate number of input arguments.
narginchk(4,5);
% Validate number of output arguments.
nargoutchk(1,1); 

% Generate Matrix D and G..................................................
% Calculate size and ValidMask.
[Ny, Nx] = size(sx);
ValidMask = isfinite(sx) & isfinite(sy);

% Expand in x-direction.
sxEx = [NaN(Ny,1), sx, NaN(Ny,2)];
xEx  = [NaN(Ny,1),  x, NaN(Ny,2)];
syEx = [NaN(Ny,1), sy, NaN(Ny,2)];
yEx  = [NaN(Ny,1),  y, NaN(Ny,2)];

ExpandMaskx = isnan(sxEx);
se = [1 1 0 1 0];
DilatedExpandMaskx = imdilate(ExpandMaskx,se);
Maskx = DilatedExpandMaskx(:,2:end-2) & ~ExpandMaskx(:,2:end-2);

% Expand in y-direction.
sxEy = [NaN(1,Nx); sx; NaN(2,Nx)];
xEy  = [NaN(1,Nx);  x; NaN(2,Nx)];
syEy = [NaN(1,Nx); sy; NaN(2,Nx)];
yEy  = [NaN(1,Nx);  y; NaN(2,Nx)];

ExpandMasky = isnan(syEy);
se = [1;1;0;1;0];
DilatedExpandMasky = imdilate(ExpandMasky,se);
Masky = DilatedExpandMasky(2:end-2,:) & ~ExpandMasky(2:end-2,:);

% Compose matrices Dx and Dy.
Num = Ny*Nx;
ee = ones(Num,1);
Dx = spdiags([-ee,ee],[0,Ny],Num,Num);
Dy = spdiags([-ee,ee],[0, 1],Num,Num);

% Compose matrices Gx and Gy.
% O(h^5)
gx_x = (-1/13*sxEx(:,1:end-3)+sxEx(:,2:end-2)+sxEx(:,3:end-1)-1/13*sxEx(:,4:end)) ...
    .*(xEx(:,3:end-1)-xEx(:,2:end-2))*13/24;
gx_y = (-1/13*syEx(:,1:end-3)+syEx(:,2:end-2)+syEx(:,3:end-1)-1/13*syEx(:,4:end)) ...
    .*(yEx(:,3:end-1)-yEx(:,2:end-2))*13/24;
gy_x = (-1/13*sxEy(1:end-3,:)+sxEy(2:end-2,:)+sxEy(3:end-1,:)-1/13*sxEy(4:end,:)) ...
    .*(xEy(3:end-1,:)-xEy(2:end-2,:))*13/24;
gy_y = (-1/13*syEy(1:end-3,:)+syEy(2:end-2,:)+syEy(3:end-1,:)-1/13*syEy(4:end,:)) ...
    .*(yEy(3:end-1,:)-yEy(2:end-2,:))*13/24;

Gx = gx_x + gx_y;
Gy = gy_x + gy_y;

% O(h^3)
gx3_x = (sxEx(:,2:end-2)+sxEx(:,3:end-1)).*(xEx(:,3:end-1)-xEx(:,2:end-2))/2;
gx3_y = (syEx(:,2:end-2)+syEx(:,3:end-1)).*(yEx(:,3:end-1)-yEx(:,2:end-2))/2;
gy3_x = (sxEy(2:end-2,:)+sxEy(3:end-1,:)).*(xEy(3:end-1,:)-xEy(2:end-2,:))/2;
gy3_y = (syEy(2:end-2,:)+syEy(3:end-1,:)).*(yEy(3:end-1,:)-yEy(2:end-2,:))/2;

Gx3 = gx3_x + gx3_y;
Gy3 = gy3_x + gy3_y;

% Use O(h^3) values, if O(h^5) is not available.
Gx(Maskx) = Gx3(Maskx);
Gy(Masky) = Gy3(Masky);

clear sx sy x y Gx3 Gy3;

% Remove NaN.
if nargin==4
    
    % Compose D.
    D = [Dx(isfinite(Gx),:); Dy(isfinite(Gy),:)];
    clear Dx Dy;
    
    % Compose G.
    G = [Gx(isfinite(Gx)); Gy(isfinite(Gy))];
    clear Gx Gy;

    % Solve "Rank deficient" for complete dataset by assuming Z(Ind)=0.  
    Ind = find(D(1,:)==-1,1);
    D(:,Ind) = [];   
    Z = D\G; 
    Z = [Z(1:Ind-1);0;Z(Ind:end)];      
    
elseif nargin==5
    
    % Compose Dz.
    Dz = spdiags(ee,0,Num,Num);
    
    % Compose D.
    D = [Dx(isfinite(Gx),:); Dy(isfinite(Gy),:); Dz(isfinite(z),:)];
    clear Dx Dy Dz;

    % Compose G.
    G = [Gx(isfinite(Gx)); Gy(isfinite(Gy)); z(isfinite(z))];
    clear Gx Gy z;
    
    % Calculate Z with least squares method.
    Z = D\G; 

end
clear D G;

% Reconstructed result.
z_hfli2q = reshape(Z,Ny,Nx);
z_hfli2q(~ValidMask)= nan;

end
