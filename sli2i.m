function z_sli2i  = sli2i(sx, sy, x, y)
%SLI2Q Spline-based Least-squares integration in quadrilateral geometry
%   (interpolation before and after).
%   D * Z = G (G is mainly composed by spline estimated values).
%
%   Reference: 
%   [1] L. Huang, J. Xue, B. Gao, C. Zuo, and M. Idir, "Spline based least 
%   squares integration for two-dimensional shape or wavefront 
%   reconstruction," Optics and Lasers in Engineering 91, 221-226 (2017)
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
zr_sli2 = sli2(sxr, syr, xr, yr);
warning('on','MATLAB:rankDeficientMatrix');

% Interpolation on height..................................................
vr = isfinite(sxr) & isfinite(syr) & isfinite(xr) & isfinite(yr);
Fzr = scatteredInterpolant(xr(vr),yr(vr),zr_sli2(vr),'natural');
z_sli2i = Fzr(x, y); 


end % End of Function.
