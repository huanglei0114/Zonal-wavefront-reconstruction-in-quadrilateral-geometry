% Integrate Slopes in Quadrilateral Geometry.

%% Clear everything.

close all
clear
clc

%% Load analytical height and slopes.

load Step_01_GenerateSlopesInQuadrilateralGeometry.mat ...
    x y sx_anal sy_anal z_anal

% Set the input slopes as analytical values.
sx = sx_anal;
sy = sy_anal;

%% Reconstruct height from slopes.

z_tfli2 = tfli2(sx, sy, x, y);
z_tfli2q = tfli2q(sx, sy, x, y);
z_hfli2 = hfli2(sx, sy, x, y);
z_hfli2q = hfli2q(sx, sy, x, y);
z_sli2 = sli2(sx, sy, x, y);
z_sli2q = sli2q(sx, sy, x, y);

t_tfli2 = timeit(@()tfli2(sx, sy, x, y));
t_tfli2q = timeit(@()tfli2q(sx, sy, x, y));
t_hfli2 = timeit(@()hfli2(sx, sy, x, y));
t_hfli2q = timeit(@()hfli2q(sx, sy, x, y));
t_sli2 = timeit(@()sli2(sx, sy, x, y));
t_sli2q = timeit(@()sli2q(sx, sy, x, y));

time = [t_tfli2 t_tfli2q; t_hfli2 t_hfli2q; t_sli2 t_sli2q];
disp('time = ');
disp(time);

[ez_tfli2, std_ez_tfli2] = EvaluateError(z_tfli2, z_anal);
[ez_tfli2q, std_ez_tfli2q] = EvaluateError(z_tfli2q, z_anal);

[ez_hfli2, std_ez_hfli2] = EvaluateError(z_hfli2, z_anal);
[ez_hfli2q, std_ez_hfli2q] = EvaluateError(z_hfli2q, z_anal);

[ez_sli2, std_ez_sli2] = EvaluateError(z_sli2, z_anal);
[ez_sli2q, std_ez_sli2q] = EvaluateError(z_sli2q, z_anal);


ez = [ez_tfli2 ez_tfli2q; ez_hfli2 ez_hfli2q; ez_sli2 ez_sli2q];
std_ez = [std_ez_tfli2, std_ez_tfli2q; std_ez_hfli2, std_ez_hfli2q; std_ez_sli2, std_ez_sli2q];
disp('std_ez = ');
disp(std_ez);

ezq = [ez_tfli2q; ez_hfli2q; ez_sli2q];
std_ezq = [std_ez_tfli2q; std_ez_hfli2q; std_ez_sli2q];
disp('std_ezq = ');
disp(std_ezq);

figure; 
subplot(121);
imshow(ez,[]); 
colorbar;
subplot(122);
imshow(ezq,[]); 
colorbar;
colormap jet;

figure; 
subplot(211);
bar(std_ez); 
subplot(212);
bar(time); 

