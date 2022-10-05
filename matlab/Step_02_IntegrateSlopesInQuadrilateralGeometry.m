% Integrate Slopes in Quadrilateral Geometry.

%% Clear everything.

close all
clear
clc

%% Load analytical height and slopes.

load Step_01_GenerateSlopesInQuadrilateralGeometry.mat ...
    x y sx_anal sy_anal z_anal

% Add normally distributed random errors to analytical slopes.
STD_slope_error = 0e-3;
sx = sx_anal + randn(size(sx_anal))*STD_slope_error;
sy = sy_anal + randn(size(sy_anal))*STD_slope_error;


%% Reconstruct height from slopes.

z_tfli2 = tfli2(sx, sy, x, y);
z_tfli2i = tfli2i(sx, sy, x, y);
z_tfli2q = tfli2q(sx, sy, x, y);

z_hfli2 = hfli2(sx, sy, x, y);
z_hfli2i = hfli2i(sx, sy, x, y);
z_hfli2q = hfli2q(sx, sy, x, y);

z_sli2 = sli2(sx, sy, x, y);
z_sli2i = sli2i(sx, sy, x, y);
z_sli2q = sli2q(sx, sy, x, y);

% Checking the error for reconstruction.
[ez_tfli2, rms_ez_tfli2] = EvaluateError(z_tfli2, z_anal);
[ez_tfli2i, rms_ez_tfli2i] = EvaluateError(z_tfli2i, z_anal);
[ez_tfli2q, rms_ez_tfli2q] = EvaluateError(z_tfli2q, z_anal);

[ez_hfli2, rms_ez_hfli2] = EvaluateError(z_hfli2, z_anal);
[ez_hfli2i, rms_ez_hfli2i] = EvaluateError(z_hfli2i, z_anal);
[ez_hfli2q, rms_ez_hfli2q] = EvaluateError(z_hfli2q, z_anal);

[ez_sli2, rms_ez_sli2] = EvaluateError(z_sli2, z_anal);
[ez_sli2i, rms_ez_sli2i] = EvaluateError(z_sli2i, z_anal);
[ez_sli2q, rms_ez_sli2q] = EvaluateError(z_sli2q, z_anal);

ez = [ez_tfli2, ez_tfli2i, ez_tfli2q;
    ez_hfli2, ez_hfli2i, ez_hfli2q;
    ez_sli2, ez_sli2i, ez_sli2q];

rms_ez = [rms_ez_tfli2, rms_ez_tfli2i, rms_ez_tfli2q;
    rms_ez_hfli2, rms_ez_hfli2i, rms_ez_hfli2q;
    rms_ez_sli2, rms_ez_sli2i, rms_ez_sli2q];

disp('rms_ez = ');
disp(rms_ez);


%% Figures for reconstruction error.

FontSize = 20;

%..........................................................................

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
saveas2('ez_tfli2.tif',600);

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
saveas2('ez_hfli2.tif',600);

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
saveas2('ez_sli2.tif',600);

%..........................................................................

figure('Name','ez_tfli2i'); 
surf(x,y,ez_tfli2i);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('error_{TFLIi} [mm]');
title('error_{TFLIi} [mm]');
view([0,90]);
axis image;
box on;
shading interp;
set(gca,'FontSize',FontSize,'FontWeight','Bold')
colormap(hsv(1024));
colorbar;
saveas2('ez_tfli2i.tif',600);

figure('Name','ez_hfli2i'); 
surf(x,y,ez_hfli2i);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('error_{HFLIi} [mm]');
title('error_{HFLIi} [mm]');
view([0,90]);
axis image;
box on;
shading interp;
set(gca,'FontSize',FontSize,'FontWeight','Bold')
colormap(hsv(1024));
colorbar;
saveas2('ez_hfli2i.tif',600);

figure('Name','ez_sli2i'); 
surf(x,y,ez_sli2i);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('error_{SLIi} [mm]');
title('error_{SLIi} [mm]');
view([0,90]);
axis image;
box on;
shading interp;
set(gca,'FontSize',FontSize,'FontWeight','Bold')
colormap(hsv(1024));
colorbar;
saveas2('ez_sli2i.tif',600);

%..........................................................................

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
saveas2('ez_tfli2q.tif',600);

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
saveas2('ez_hfli2q.tif',600);

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
saveas2('ez_sli2q.tif',600);


%% Figures for RMS of reconstruction error

FontSize = 16;
figure('Name','RMS of reconstruction error' ...
    ,'position',[100 200 1200 500]...
    ); 
bar(rms_ez');
colormap(lines(3));
hold on;
method_str = {'TFLI','TFLIi','TFLIq';'HFLI','HFLIi','HFLIq';'SLI','SLIi','SLIq'};
for num = 1 : numel(method_str)
    str = num2str(rms_ez(num),'%.e');
    str = [strrep(str, 'e-0', '×10^{-') '}'];
    [row,col] = ind2sub([3,3],num);
    text(col+0.22*(row-2),rms_ez(num), str ...
        ,'VerticalAlignment','Bottom' ...
        ,'HorizontalAlignment','Center'...
        ,'FontSize', FontSize ...
        ,'FontWeight', 'Bold' ...
        ,'Color', 'r' ...
        );
    text(col+0.22*(row-2),0, method_str{row,col} ...
        ,'VerticalAlignment','Top' ...
        ,'HorizontalAlignment','Center'...
        ,'FontSize', FontSize ...
        ,'FontWeight', 'Bold' ...
        );
end
box off;
set(gca,'XTick', [] ...
    ,'FontSize', FontSize ...
    ,'FontWeight', 'Bold' ...
    );
ylabel('RMSE [mm]')
saveas2('RMS of reconstruction error.tif',600);


%% Figures for computing time

if verLessThan('matlab','R2013b')
    % -- Code to run in MATLAB R2013b and earlier here --
    return
else
    % -- Code to run in MATLAB R2013b and later here --
    
    % Checking the time for reconstruction.
    t_tfli2 = timeit(@()tfli2(sx, sy, x, y));
    t_tfli2i = timeit(@()tfli2i(sx, sy, x, y));
    t_tfli2q = timeit(@()tfli2q(sx, sy, x, y));
    
    t_hfli2 = timeit(@()hfli2(sx, sy, x, y));
    t_hfli2i = timeit(@()hfli2i(sx, sy, x, y));
    t_hfli2q = timeit(@()hfli2q(sx, sy, x, y));
    
    t_sli2 = timeit(@()sli2(sx, sy, x, y));
    t_sli2i = timeit(@()sli2i(sx, sy, x, y));
    t_sli2q = timeit(@()sli2q(sx, sy, x, y));

    computing_time = [t_tfli2, t_tfli2i, t_tfli2q;
        t_hfli2, t_hfli2i, t_hfli2q;
        t_sli2, t_sli2i, t_sli2q];
    
    disp('computing_time = ');
    disp(computing_time);

    % Plot the figures.....................................................
    figure('Name','Time' ...
        ,'position',[100 200 1100 500]...
        ); 
    bar(computing_time'); 
    method_str = {'TFLI','HFLI','SLI';
        'TFLIi','HFLIi','SLIi';
        'TFLIq','HFLIq','SLIq'};
    ct = computing_time;
    colormap(lines(3));
    for num = 1 : numel(method_str)
        str = num2str(ct(num),'%.3g');
        [row,col] = ind2sub([3,3],num);
        text(col+0.11*((row-2)*2),ct(num), str ...
            ,'VerticalAlignment','Bottom' ...
            ,'HorizontalAlignment','Center'...
            ,'FontSize', FontSize ...
            ,'FontWeight', 'Bold' ...
            );
        text(col+0.11*((row-2)*2),0, method_str{col,row} ...
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
    box off;
    ylabel('Time [s]')
end

