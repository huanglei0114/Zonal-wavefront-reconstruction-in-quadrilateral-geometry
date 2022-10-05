function z_hfli2i  = hfli2i(sx, sy, x, y)
%HFLI2I Higher-order Finite-difference-based Least-squares Integration in 
%   quadrilateral geometry (interpolation before and after).
%   D * Z = G.
%
%   Reference: 
%   [1] Guanghui Li, Yanqiu Li, Ke Liu, Xu Ma, and Hai Wang, "Improving 
%   wavefront reconstruction accuracy by using integration equations with 
%   higher-order truncation errors in the Southwell geometry," J. Opt. Soc.
%   Am. A 30, 1448-1459 (2013)
%
%   [2] H. Ren, F. Gao, and X. Jiang, "Improvement of high-order 
%   least-squares integration method for stereo deflectometry," Appl. Opt. 
%   54, 10249-10255 (2015).

%   Copyright since 2016 by Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%   2016-10-01 Original Version

% Check the number of arguments............................................
% Validate number of input arguments.
narginchk(4,4);
% Validate number of output arguments.
nargoutchk(1,1); 


% Rectangular mesh for interpolation.......................................
v = isfinite(sx) & isfinite(sy) & isfinite(x) & isfinite(y);
sz = size(x);
[xr, yr] = meshgrid(linspace(min(x(v)), max(x(v)), sz(2)) ...
    , linspace(min(y(v)), max(y(v)), sz(1)));


% Interpolation on slopes..................................................
Fsx = scatteredInterpolant(x(v), y(v), sx(v), 'natural', 'none');
sxr = Fsx(xr, yr); 

Fsy = scatteredInterpolant(x(v), y(v), sy(v), 'natural', 'none');
syr = Fsy(xr, yr); 


% Integration..............................................................
warning('off','MATLAB:rankDeficientMatrix');
zr_hfli2 = hfli2(sxr, syr, xr, yr);
warning('on','MATLAB:rankDeficientMatrix');


% Interpolation on height..................................................
vr = isfinite(sxr) & isfinite(syr) & isfinite(xr) & isfinite(yr);
Fzr = scatteredInterpolant(xr(vr), yr(vr), zr_hfli2(vr), 'natural');
z_hfli2i = Fzr(x, y); 


end % End of Function.
