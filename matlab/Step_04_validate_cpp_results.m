clear;
close all;
clc;


%% load C++ results
X = Load2DMatrix('../data/X.bin', 'double');
Y = Load2DMatrix('../data/Y.bin', 'double');
Z_calc = Load2DMatrix('../data/Z_calc.bin', 'double');
Z = Load2DMatrix('../data/Z.bin', 'double');

figure;
subplot(1, 3, 1);
show_surface_map(X, Y, Z, 0, 'jet', 'flat', true, 2, 1, '', 'Nominal height');
subplot(1, 3, 2);
show_surface_map(X, Y, Z_calc, 0, 'jet', 'flat', true, 2, 1, '', 'Integrated height');
subplot(1, 3, 3);
show_surface_map(X, Y, Z - Z_calc, 0, 'jet', 'flat', true, 5, 1, '', 'Residual height');


g = Load2DMatrix('../data/g.bin', 'double');
