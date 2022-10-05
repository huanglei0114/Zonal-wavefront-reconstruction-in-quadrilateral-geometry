function z_hfli2 = hfli2(sx, sy, x, y, z)
%HFLI2 Higher-order Finite-difference-based Least-squares Integration. 
%   (Algorithm 1 in Reference)
%   D * Z = G.
%   
%   Reference: 
%   Guanghui Li, Yanqiu Li, Ke Liu, Xu Ma, and Hai Wang, "Improving 
%   wavefront reconstruction accuracy by using integration equations with 
%   higher-order truncation errors in the Southwell geometry," J. Opt. Soc.
%   Am. A 30, 1448-1459 (2013) 

%   Copyright since 2013 by Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%   2013-12-01 Original Version
%   2014-01-31 First revision on matrix
%   2014-11-13 revision on speed
%   2016-01-29 Solve "Warning: Rank deficient" for complete dataset.

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
sx = [NaN(Ny,1),sx,NaN(Ny,2)];
x  = [NaN(Ny,1),x ,NaN(Ny,2)];

ExpandMaskx = isnan(sx);
se = [1 1 0 1 0];
DilatedExpandMaskx = imdilate(ExpandMaskx,se);
Maskx = DilatedExpandMaskx(:,2:end-2) & ~ExpandMaskx(:,2:end-2);

% Expand in y-direction.
sy = [NaN(1,Nx);sy;NaN(2,Nx)];
y  = [NaN(1,Nx);y ;NaN(2,Nx)];

ExpandMasky = isnan(sy);
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
Gx = (-1/13*sx(:,1:end-3)+sx(:,2:end-2)+sx(:,3:end-1)-1/13*sx(:,4:end)) ...
    .*(x(:,3:end-1)-x(:,2:end-2))*13/24;
Gy = (-1/13*sy(1:end-3,:)+sy(2:end-2,:)+sy(3:end-1,:)-1/13*sy(4:end,:)) ...
    .*(y(3:end-1,:)-y(2:end-2,:))*13/24;

% O(h^3)
Gx3 = (sx(:,2:end-2)+sx(:,3:end-1)).*(x(:,3:end-1)-x(:,2:end-2))/2;
Gy3 = (sy(2:end-2,:)+sy(3:end-1,:)).*(y(3:end-1,:)-y(2:end-2,:))/2;

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
z_hfli2 = reshape(Z,Ny,Nx);
z_hfli2(~ValidMask)= nan;

end