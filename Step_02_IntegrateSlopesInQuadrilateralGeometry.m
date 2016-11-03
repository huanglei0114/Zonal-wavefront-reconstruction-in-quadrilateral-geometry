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

% Checking the error for reconstruction.
[ez_tfli2, rms_ez_tfli2] = EvaluateError(z_tfli2, z_anal);
[ez_tfli2q, rms_ez_tfli2q] = EvaluateError(z_tfli2q, z_anal);

[ez_hfli2, rms_ez_hfli2] = EvaluateError(z_hfli2, z_anal);
[ez_hfli2q, rms_ez_hfli2q] = EvaluateError(z_hfli2q, z_anal);

[ez_sli2, rms_ez_sli2] = EvaluateError(z_sli2, z_anal);
[ez_sli2q, rms_ez_sli2q] = EvaluateError(z_sli2q, z_anal);

ez = [ez_tfli2 ez_tfli2q; ez_hfli2 ez_hfli2q; ez_sli2 ez_sli2q];
rms_ez = [rms_ez_tfli2, rms_ez_tfli2q; rms_ez_hfli2, rms_ez_hfli2q; rms_ez_sli2, rms_ez_sli2q];
disp('rms_ez = ');
disp(rms_ez);

%% Figures for reconstruction error.
FontSize = 20;

figure('Name','ez_tfli2'); 
surf(x,y,ez_tfli2);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('error_{TFLI} [mm]');
title('error_{TFLI} [mm]');
view([0,90]);
axis image;
box on;
shading interp;
set(gca,'FontSize',FontSize,'FontWeight','Bold')
colormap(hsv(1024));
colorbar;

figure('Name','ez_hfli2'); 
surf(x,y,ez_hfli2);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('error_{HFLI} [mm]');
title('error_{HFLI} [mm]');
view([0,90]);
axis image;
box on;
shading interp;
set(gca,'FontSize',FontSize,'FontWeight','Bold')
colormap(hsv(1024));
colorbar;

figure('Name','ez_sli2'); 
surf(x,y,ez_sli2);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('error_{SLI} [mm]');
title('error_{SLI} [mm]');
view([0,90]);
axis image;
box on;
shading interp;
set(gca,'FontSize',FontSize,'FontWeight','Bold')
colormap(hsv(1024));
colorbar;

figure('Name','ez_tfli2q'); 
surf(x,y,ez_tfli2q);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('error_{TFLIq} [mm]');
title('error_{TFLIq} [mm]');
view([0,90]);
axis image;
box on;
shading interp;
set(gca,'FontSize',FontSize,'FontWeight','Bold')
colormap(hsv(1024));
colorbar;
zlim([-1,1])

figure('Name','ez_hfli2q'); 
surf(x,y,ez_hfli2q);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('error_{HFLIq} [mm]');
title('error_{HFLIq} [mm]');
view([0,90]);
axis image;
box on;
shading interp;
set(gca,'FontSize',FontSize,'FontWeight','Bold')
colormap(hsv(1024));
colorbar;

figure('Name','ez_sli2q'); 
surf(x,y,ez_sli2q);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('error_{SLIq} [mm]');
title('error_{SLIq} [mm]');
view([0,90]);
axis image;
box on;
shading interp;
set(gca,'FontSize',FontSize,'FontWeight','Bold')
colormap(hsv(1024));
colorbar;


%% figures for RMS of reconstruction error

FontSize = 16;
figure('Name','RMS of reconstruction error' ...
    ,'position',[100 200 800 300]...
    ); 
bar(rms_ez');
hold on;
method_str = {'TFLI','TFLIq';'HFLI','HFLIq';'SLI','SLIq'};
for num = 1 : 6
    str = num2str(rms_ez(num),'%.e');
    str = [strrep(str, 'e-0', '×10^{-') '}'];
    [row,col] = ind2sub([3,2],num);
    text(col+0.225*(row-2),rms_ez(num), str ...
        ,'VerticalAlignment','Bottom' ...
        ,'HorizontalAlignment','Center'...
        ,'FontSize', FontSize ...
        ,'FontWeight', 'Bold' ...
        ,'Color', 'r' ...
        );
    text(col+0.225*(row-2),0, method_str{row,col} ...
        ,'VerticalAlignment','Top' ...
        ,'HorizontalAlignment','Center'...
        ,'FontSize', FontSize ...
        ,'FontWeight', 'Bold' ...
        );
end

set(gca,'XTick', [] ...
    ,'FontSize', FontSize ...
    ,'FontWeight', 'Bold' ...
    );
ylabel('RMSE [mm]')



%% figures for computing time

if verLessThan('matlab','R2013b')
    % -- Code to run in MATLAB R2013b and earlier here --
    return
else
    % -- Code to run in MATLAB R2013b and later here --
    % Checking the time for reconstruction.
    t_tfli2 = timeit(@()tfli2(sx, sy, x, y));
    t_tfli2q = timeit(@()tfli2q(sx, sy, x, y));
    t_hfli2 = timeit(@()hfli2(sx, sy, x, y));
    t_hfli2q = timeit(@()hfli2q(sx, sy, x, y));
    t_sli2 = timeit(@()sli2(sx, sy, x, y));
    t_sli2q = timeit(@()sli2q(sx, sy, x, y));

    computing_time = [t_tfli2 t_tfli2q; t_hfli2 t_hfli2q; t_sli2 t_sli2q];
    disp('computing_time = ');
    disp(computing_time);

    figure('Name','Time' ...
        ,'position',[100 200 800 300]...
        ); 
    bar(computing_time); 
    method_str = {'TFLI','HFLI','SLI';'TFLIq','HFLIq','SLIq'};
    ct = computing_time';
    color = parula(2);
    for num = 1 : 6
        str = num2str(ct(num),'%.3f');
        [row,col] = ind2sub([2,3],num);
        text(col+0.15*((row-1.5)*2),ct(num), str ...
            ,'VerticalAlignment','Bottom' ...
            ,'HorizontalAlignment','Center'...
            ,'FontSize', FontSize-4 ...
            ,'FontWeight', 'Bold' ...
            );
        text(col+0.15*((row-1.5)*2),0, method_str{row,col} ...
            ,'VerticalAlignment','Top' ...
            ,'HorizontalAlignment','Center'...
            ,'FontSize', FontSize ...
            ,'FontWeight', 'Bold' ...
            );
    end

    set(gca,'XTickLabel', []...
        ,'FontSize', FontSize ...
        ,'FontWeight', 'Bold' ...
        )

    ylabel('Time [s]')
end

