function z_sli2 = sli2(sx, sy, x, y)
%SLI2 Spline-based Least-squares integration.
%   D * Z = G (G is mainly composed by spline estimated values).

%   Copyright since 2016 by Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%   2016-09-29 Original Version
%   2016-11-01 Revised for x and y increasing directions.

% Check the number of arguments............................................
% Validate number of input arguments.
narginchk(4,4);
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
ee = ones(Ny*Nx,1);
Dx = spdiags([-ee,ee],[0,Ny],Ny*(Nx-1),Ny*Nx); 
Dy = spdiags([-ee,ee],[0,1],Ny*Nx,Ny*Nx);

% Compose matrices Gx and Gy.
Gx = (sx(:,1:end-1)+sx(:,2:end)).*(x(:,2:end)-x(:,1:end-1))/2;  
Gy = (sy(1:end-1,:)+sy(2:end,:)).*(y(2:end,:)-y(1:end-1,:))/2; 

% Compose D.
D = [Dx(isfinite(Gx),:); Dy(isfinite(Gy),:)];
clear Dx Dy;

% Compose matrix SpGx.
spGx = ComposeSpGx(x,sx,ValidMask,Nx,Ny);
% Compose matrix SpGy.
spGy = ComposeSpGy(y,sy,ValidMask,Nx,Ny);
clear sx sy x y;

% Replace with spline values, if available.
Gy(end,:)=[];
Gx(isfinite(spGx)) = spGx(isfinite(spGx));
Gy(isfinite(spGy)) = spGy(isfinite(spGy));

% Compose G.
G = [Gx(isfinite(Gx)); Gy(isfinite(Gy))];
clear Gx Gy;

% Solve "Warning: Rank deficient" for complete data by assuming Z(Ind)=0.  
Ind = find(D(1,:)==-1,1);
D(:,Ind) = [];
Z = D\G;
Z = [Z(1:Ind-1);0;Z(Ind:end)];

% Reconstructed result.
z_sli2 = reshape(Z,Ny,Nx);
z_sli2(~ValidMask)= nan;

end




%% Subfunctions.

% Compose matrix spGx
function SpGx = ComposeSpGx(x,sx,ValidMask,Nx,Ny)

SpGx = NaN(Ny,Nx-1);
for ny = 1:Ny
    xl = x(ny,:)';
    vl = sx(ny,:)';
    
    % Check the number of sections.
    ml = ValidMask(ny,:)';   
    [Ns, Indices] = CheckSection(ml, Nx);
    
    % Spline fitting section by section.
    gs = cell(Ns,1);
    for ns = 1:Ns
        idx = Indices{ns};
        xx = xl(idx);
        vv = vl(idx);
        if length(xx)>1
            pp = spline(xx,vv); % "not-a-knot end condition"
            c = pp.coefs;
            switch(size(c,2))
                case 4  % 4 points for piecewise cubic spline fitting.
                    dx = diff(xx);
                    if sign(mean(dx))==1
                        gs{ns} = dx.*(c(:,4) + dx.*(c(:,3)./2 + dx.*(c(:,2)./3 + dx.*c(:,1)./4)));
                    else
                        dx = -flipud(dx);
                        gs{ns} = dx.*(c(:,4) + dx.*(c(:,3)./2 + dx.*(c(:,2)./3 + dx.*c(:,1)./4)));
                        gs{ns} = -flipud(gs{ns});
                    end 
                    
                case 3  % 3 points for 2nd order polynominal fitting.
                    % Here we do not use the polynomials. 
                    % We are going to use the Southwell expression instead
                    gs{ns} = diff(xx)*NaN; 
                
                case 2  % 2 points for 1st order polynominal fitting.
                    % Here we do not use the polynomials. 
                    % We are going to use the Southwell expression instead
                    gs{ns} = diff(xx)*NaN;
                
                case 1
                    % Logically impossible.
                    error('Only one point for fitting in x direction!');
                
                otherwise
                    % Logically impossible.
                    error('Unexpected number of points for fitting in x direction!');
            end
        end
    end
    sg = cat(1,gs{:});
    Valid = ml(1:end-1) & ml(1+1:end);
    pt = 1;
    for nx = 1 : Nx-1
        if Valid(nx) == 1
            SpGx(ny, nx) = sg(pt);
            pt = pt + 1;
        end
    end
end
end


% Compose matrix spGy
function SpGy = ComposeSpGy(y,sy,ValidMask,Nx,Ny)
SpGy = NaN(Ny-1,Nx);
for nx = 1:Nx
    yl = y(:,nx);
    vl = sy(:,nx);
    
    % Check the number of sections.
    ml = ValidMask(:,nx);
    [Ns, Indices] = CheckSection(ml, Ny);
    
    % Spline fitting section by section.
    gs = cell(Ns,1);
    for ns = 1:Ns
        idx = Indices{ns};
        yy = yl(idx);
        vv = vl(idx);
        if length(yy)>1
            pp = spline(yy,vv); % "not-a-knot end condition"
            c = pp.coefs;
            switch(size(c,2))
                case 4  % 4 points for piecewise cubic spline fitting.
                    dy = diff(yy);
                    if sign(mean(dy))==1
                        gs{ns} = dy.*(c(:,4) + dy.*(c(:,3)./2 + dy.*(c(:,2)./3 + dy.*c(:,1)./4)));
                    else
                        dy = -flipud(diff(yy));
                        gs{ns} = dy.*(c(:,4) + dy.*(c(:,3)./2 + dy.*(c(:,2)./3 + dy.*c(:,1)./4)));
                        gs{ns} = -flipud(gs{ns});
                    end
                    
                case 3  % 3 points
                    % Here we do not use the polynomials. 
                    % We are going to use the Southwell expression instead
                    gs{ns} = diff(yy)*NaN;

                case 2  % 2 points
                    % Here we do not use the polynomials. 
                    % We are going to use the Southwell expression instead
                    gs{ns} = diff(yy)*NaN;
                
                case 1
                    % Logically impossible.
                    error('Only one point for fitting in y direction!');
                
                otherwise
                    % Logically impossible.
                    error('Unexpected number of points for fitting in y direction!');
            end
        end
    end
    sg = cat(1,gs{:});
    Valid = ml(1:end-1) & ml(1+1:end);
    pt = 1;
    for ny = 1 : Ny-1
        if Valid(ny) == 1
            SpGy(ny, nx) = sg(pt);
            pt = pt + 1;
        end
    end
end
end


% Check Sections.
function [Ns, Indices] = CheckSection(ml, N)
if all(ml)==true      
    Ns = 1;
    Indices{Ns} = 1:N;
else
    Indices = cell(N,1);
    first = nan;
    last = nan;
    Ns = 0;
    for n = 1:N
        % Find the first.
        if n==1
            if ml(n)==true
                first = n;
            end
        else
            if ml(n)==true && ml(n-1)==false
                first = n;
            end
        end

        % Find the last.
        if n==N
            if ml(n)==true
                last = n;
            end
        else
            if ml(n)==true && ml(n+1)==false
                last = n;
            end
        end

        % Sum up the total number of sections and compose the Indices.
        if isfinite(first) && isfinite(last)
            Ns = Ns + 1;
            Indices{Ns} = first:last;
            first = nan;
            last = nan;
        end
    end
end
end
