clear;
close all;
clc;


%% load C++ results
X = Load2DMatrix('../data/X.bin', 'double');
Y = Load2DMatrix('../data/Y.bin', 'double');
Z_calc = Load2DMatrix('../data/Z_calc.bin', 'double');

show_surface_map(X, Y, Z_calc, 0, 'jet', 'flat', true, 2, 1, '', '');

g = Load2DMatrix('../data/g.bin', 'double');
