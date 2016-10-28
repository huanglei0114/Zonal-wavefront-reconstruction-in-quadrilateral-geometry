function z_tfli2  = tfli2(sx, sy, x, y, z)
%TFLI2 Traditional Finite-difference-based Least-squares Integration.
%   D * Z = G.
% 
%   Reference:
%   W.H. Southwell, "Wave-front estimation from wave-front slope 
%   measurements," J. Opt. Soc. Am. 70, 998-1006 (1980) 

%   Copyright since 2010 by Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%   2010-10-01 Original Version
%   2014-02-01 Revision on data structure

% Check the number of arguments............................................
% Validate number of input arguments.
narginchk(4,5);
% Validate number of output arguments.
nargoutchk(1,1); 

% Generate Matrix D and G..................................................
% Calculate size and ValidMask.
[Ny, Nx] = size(sx);
ValidMask = isfinite(sx) & isfinite(sy);

% Expand sy and y.
sy = [sy;NaN(1,Nx)];
y  = [y ;NaN(1,Nx)];

% Compose matrices Dx and Dy.
Num = Ny*Nx;
ee = ones(Num,1);
Dx = spdiags([-ee,ee],[0,Ny],Ny*(Nx-1),Num); 
Dy = spdiags([-ee,ee],[0, 1],Num,Num);

% Compose matrices Gx and Gy.
Gx = (sx(:,1:end-1)+sx(:,2:end)).*(x(:,2:end)-x(:,1:end-1))/2;  
Gy = (sy(1:end-1,:)+sy(2:end,:)).*(y(2:end,:)-y(1:end-1,:))/2; 

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
z_tfli2 = reshape(Z,Ny,Nx);
z_tfli2(~ValidMask)= nan;

end % End of Function.

