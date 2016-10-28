function z_tfli2q  = tfli2q(sx, sy, x, y, z)
%TFLI2 Traditional Finite-difference-based Least-squares Integration in 
%   quadrilateral geometry.
%   D * Z = G.
% 
%   Reference:
%   W.H. Southwell, "Wave-front estimation from wave-front slope 
%   measurements," J. Opt. Soc. Am. 70, 998-1006 (1980) 

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

% Expand in y direction.
sxEy = [sx; NaN(1,Nx)];
xEy  = [ x; NaN(1,Nx)];
syEy = [sy; NaN(1,Nx)];
yEy  = [ y; NaN(1,Nx)];

% Compose matrices Dx and Dy.
Num = Ny*Nx;
ee = ones(Num,1);
Dx = spdiags([-ee,ee],[0,Ny],Ny*(Nx-1),Num); 
Dy = spdiags([-ee,ee],[0, 1],Num,Num);

% Compose matrices Gx and Gy.
gx_x = (sx(:,1:end-1)+sx(:,2:end)).*(x(:,2:end)-x(:,1:end-1))/2;
gx_y = (sy(:,1:end-1)+sy(:,2:end)).*(y(:,2:end)-y(:,1:end-1))/2;  
gy_x = (sxEy(1:end-1,:)+sxEy(2:end,:)).*(xEy(2:end,:)-xEy(1:end-1,:))/2;
gy_y = (syEy(1:end-1,:)+syEy(2:end,:)).*(yEy(2:end,:)-yEy(1:end-1,:))/2;  
Gx = gx_x + gx_y;
Gy = gy_x + gy_y;

clear sx sy x y;

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
z_tfli2q = reshape(Z,Ny,Nx);
z_tfli2q(~ValidMask)= nan;

end % End of Function.

