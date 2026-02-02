BASE_DIRECTORY = "../..";

addpath(BASE_DIRECTORY);

% ================================================================
% PARAMS
% ================================================================

HISTOGRAM_BINCOUNT = 180;
HISTOGRAM_CLIM = [0.2 0.6] * 10^-3;

RADIAL_HISTOGRAM_BINCOUNT = 80;
RADIAL_HISTOGRAM_USEPLIM = false;
RADIAL_HISTOGRAM_PLIM = 0.025;

TRAJECTORY_CLIM = [25, 45];

SPEED_HISTOGRAM_BINCOUNT = RADIAL_HISTOGRAM_BINCOUNT;
SPEED_HISTOGRAM_MAXSPD = 40;

WAVEFIELD_CLIM = [-1 +1] * 0.010;

CROSS_SECTIONS_DROPLET_SIZE = 20;

SPATIALSPEED_RES_BASE = 180;
SPATIALSPEED_BINSCALE = 1;
SPATIALSPEED_INTERP_SCALE = 4;
SPATIALSPEED_CLIM = [14, 25];

% Saving
VAR_outputFolder = "RES";

%% ================================================================

folders = dir(fullfile(VAR_outputFolder, "RES_*"));
threadCount = length(folders);

ISLOCAL = isfile(BASE_DIRECTORY + "/ISLOCAL");

% Only do one run if using on local
if ISLOCAL
  threadCount = 1;
end

parfor i = 1:threadCount
    folderPath = fullfile(folders(i).folder, folders(i).name);

    % pathId is important for saving!
    [pathstr, name, ext] = fileparts(folderPath);
    pathId = fullfile(folderPath, strcat(name, ext)); % Add duplicate of innermost folder string

    fileSearchPath = fullfile(folderPath, "*.mat");
    files = dir(fileSearchPath);

    x_data = [];
    y_data = [];
    x_runs = cell(numel(files), 1);
    y_runs = cell(numel(files), 1);
    for j = 1:numel(files)
        filePath = fullfile(files(j).folder, files(j).name);
        loadP = load(filePath);
        p = loadP.p;

        x_data = [x_data; p.x_data];
        y_data = [y_data; p.y_data];

        x_runs{j} = p.x_data;
        y_runs{j} = p.y_data;
    end

    close all

    %% Setup
    xs = [];
    ys = [];
    for j = 1:size(x_data, 2)
        xs = [xs; x_data(:, j)];
        ys = [ys; y_data(:, j)];
    end

    BOUNDS_R_PERL = [-0.5 0.5] * p.Lx;
    BOUNDS_R_PERL_CROPPED = [-(p.Rc + 1), (p.Rc + 1)];
    
    %% HISTOGRAM
    figure
    
    xEdges = linspace(-p.Lx/2, p.Lx/2, HISTOGRAM_BINCOUNT + 1);
    yEdges = linspace(-p.Ly/2, p.Ly/2, HISTOGRAM_BINCOUNT + 1);
    [bins, xEdges, yEdges] = histcounts2(xs, ys, xEdges, yEdges, Normalization="probability");

    axisValuesX = linspace(-p.Lx/2, p.Lx/2, HISTOGRAM_BINCOUNT);
    axisValuesY = linspace(-p.Ly/2, p.Ly/2, HISTOGRAM_BINCOUNT);
    hold on
    imagesc(axisValuesX, axisValuesY, bins');
    DrawBounds2D(p)
    hold off
    
    title("Histogram", 'Interpreter', 'latex')
    xlabel('$x/\lambda_F$', Interpreter = "latex")
    ylabel('$y/\lambda_F$', Interpreter = "latex")
    axis square
    colorbar
    xlim(BOUNDS_R_PERL_CROPPED)
    ylim(BOUNDS_R_PERL_CROPPED)
    clim(HISTOGRAM_CLIM)

    exportgraphics(gca, pathId + "_histogram.png");
    
    %% RADIAL HISTOGRAM
    figure
    
    r_data = sqrt(xs.^2 + ys.^2);
    radEdges = linspace(0, p.Rc, RADIAL_HISTOGRAM_BINCOUNT);
    radBins = histcounts(r_data, radEdges);
    radCenters = (radEdges(1:end - 1) + radEdges(2:end)) / 2;
    radProb = radBins ./ radCenters;
    radProb = radProb / length(r_data);
    
    j0Domain = (0:0.02:p.Rc);
    J0 = besselj(0, (2*pi) * j0Domain);
    
    bar(radCenters, radProb, 'hist');
    hold on
    plot(j0Domain, abs(J0) * max(radProb), Color="red");
    DrawBounds1DRad(p)
    hold off

    title("Radial Histogram", 'Interpreter', 'latex')
    xlabel('$r/\lambda_F$', Interpreter = "latex")
    ylabel('$p$', Interpreter = "latex")
    xlim([0 p.Rc]);
    if (RADIAL_HISTOGRAM_USEPLIM)
        ylim([0 RADIAL_HISTOGRAM_PLIM]);
    end

    exportgraphics(gca, pathId + "_radial_histogram.png");

    %% TRAJECTORY
    figure
    axes = gca;
    hold on

    for j = 1:size(p.x_data, 2)
        xs_pl = p.x_data(:, j);
        ys_pl = p.y_data(:, j);
    
        vx_plpf = diff(xs_pl, 1, 1);
        vx_plpf = mod(vx_plpf + p.Lx/2, p.Lx) - p.Lx/2;
    
        vy_plpf = diff(ys_pl, 1, 1);
        vy_plpf = mod(vy_plpf + p.Ly/2, p.Ly) - p.Ly/2;
    
        vs_mmpf = sqrt(vx_plpf.^2 + vy_plpf.^2) * p.lambdaF * 1000;
        vs_mmps = vs_mmpf * 1/p.TF;

        % Need to drop one so xs.Count and ys.Count match vs.Count
        xs_pl = xs_pl(2:end, :);
        ys_pl = ys_pl(2:end, :);

        cline(xs_pl, ys_pl, vs_mmps);
    end
    DrawBounds2D(p)

    title("Trajectories", 'Interpreter', 'latex')
    xlabel('$x/\lambda_F$', Interpreter = "latex")
    ylabel('$y/\lambda_F$', Interpreter = "latex")
    axis square
    xlim(BOUNDS_R_PERL_CROPPED)
    ylim(BOUNDS_R_PERL_CROPPED)
    clim(TRAJECTORY_CLIM)

    exportgraphics(axes, pathId + "_trajectories.png");
    
    %% SPEED HISTOGRAM
    figure

    velEdges = linspace(0, SPEED_HISTOGRAM_MAXSPD, SPEED_HISTOGRAM_BINCOUNT + 1);
    velBins = zeros(1, SPEED_HISTOGRAM_BINCOUNT);
    for j = 1:numel(x_runs)
        xs_pl = x_runs{j};
        ys_pl = y_runs{j};
    
        vx_plpf = diff(xs_pl, 1, 1);
        vx_plpf = mod(vx_plpf + p.Lx/2, p.Lx) - p.Lx/2;
    
        vy_plpf = diff(ys_pl, 1, 1);
        vy_plpf = mod(vy_plpf + p.Ly/2, p.Ly) - p.Ly/2;
    
        vs_mmpf = sqrt(vx_plpf.^2 + vy_plpf.^2) * p.lambdaF * 1000;
        vs_mmps = vs_mmpf * 1/p.TF;
    
        counts = histcounts(vs_mmps, velEdges);
    
        velBins = velBins + counts;
    end

    histogram(BinEdges = velEdges, BinCounts = velBins)

    title("Speed Histogram", 'Interpreter', 'latex')
    xlabel('$v$ (mm/s)', Interpreter = "latex")
    ylabel('$p$', Interpreter = "latex")
    xlim([0 SPEED_HISTOGRAM_MAXSPD])

    exportgraphics(gca, pathId + "_speed_histogram.png");
    
    %% SPEED PLOT
    figure

    xs_pl = p.x_data;
    ys_pl = p.y_data;

    vx_plpf = diff(xs_pl, 1, 1);
    vx_plpf = mod(vx_plpf + p.Lx/2, p.Lx) - p.Lx/2;

    vy_plpf = diff(ys_pl, 1, 1);
    vy_plpf = mod(vy_plpf + p.Ly/2, p.Ly) - p.Ly/2;

    vs_mmpf = sqrt(vx_plpf.^2 + vy_plpf.^2) * p.lambdaF * 1000;
    vs_mmps = vs_mmpf * 1/p.TF;

    plot(vs_mmps);

    title("Speed", 'Interpreter', 'latex')
    xlabel('$t_n$', Interpreter = "latex")
    ylabel('$v$ (mm/s)', Interpreter = "latex")

    exportgraphics(gca, pathId + "_speed_plot.png");
    
    %% SPATIALSPEED
    figure

    SPATIALSPEED_RES = SPATIALSPEED_RES_BASE * SPATIALSPEED_BINSCALE;
    SPATIALSPEED_BINCOUNT = SPATIALSPEED_RES;

    xEdges = linspace(-p.Lx/2, p.Lx/2, SPATIALSPEED_RES + 1);
    yEdges = linspace(-p.Ly/2, p.Ly/2, SPATIALSPEED_RES + 1);
    axisValuesX = linspace(-p.Lx/2, p.Lx/2, SPATIALSPEED_RES);
    axisValuesY = linspace(-p.Ly/2, p.Ly/2, SPATIALSPEED_RES);

    xs_pl = p.x_data;
    ys_pl = p.y_data;

    vx_plpf = diff(xs_pl, 1, 1);
    vx_plpf = mod(vx_plpf + p.Lx/2, p.Lx) - p.Lx/2;

    vy_plpf = diff(ys_pl, 1, 1);
    vy_plpf = mod(vy_plpf + p.Ly/2, p.Ly) - p.Ly/2;

    vs_mmpf = sqrt(vx_plpf.^2 + vy_plpf.^2) * p.lambdaF * 1000;
    vs_mmps = vs_mmpf * 1/p.TF;

    % Need to drop one so xs.Count and ys.Count match vs.Count
    ss_xs_pl = xs_pl(2:end, :);
    ss_ys_pl = ys_pl(2:end, :);
    
    % Populate bins
    bins = nan(SPATIALSPEED_BINCOUNT, SPATIALSPEED_BINCOUNT);
    for xIdx = 1:SPATIALSPEED_BINCOUNT
        for yIdx = 1:SPATIALSPEED_BINCOUNT
            in_bin = ss_xs_pl >= xEdges(xIdx) & ss_xs_pl < xEdges(xIdx+1) & ...
                     ss_ys_pl >= yEdges(yIdx) & ss_ys_pl < yEdges(yIdx+1);
            if any(in_bin)
                bins(xIdx, yIdx) = mean(vs_mmps(in_bin));
            end
        end
    end

    hold on
    imagesc(axisValuesX, axisValuesY, bins');
    DrawBounds2D(p)
    hold off

    title("Spatial Speed", 'Interpreter', 'latex')
    xlabel('$x/\lambda_F$', Interpreter = "latex")
    ylabel('$y/\lambda_F$', Interpreter = "latex")
    axis square
    colorbar
    xlim(BOUNDS_R_PERL_CROPPED)
    ylim(BOUNDS_R_PERL_CROPPED)
    clim(SPATIALSPEED_CLIM)

    exportgraphics(gca, pathId + "_spatialspeed.png");

    %% WAVEFIELD
    figure
    wavefield = p.eta_data(:, :, end);
    DrawWavefield2D(p, wavefield);
    hold on
    viscircles([p.x_data(end,:); p.y_data(end,:)]', p.drop_radius / p.lambdaF * ones(1, p.n_drops));
    hold off
    
    title("Wavefield", 'Interpreter', 'latex')
    xlabel('$x/\lambda_F$', Interpreter = "latex")
    ylabel('$y/\lambda_F$', Interpreter = "latex")
    axis square
    xlim(BOUNDS_R_PERL_CROPPED)
    ylim(BOUNDS_R_PERL_CROPPED)
    colorbar
    clim(WAVEFIELD_CLIM)

    exportgraphics(gca, pathId + "_wavefield.png");

    %% AVERAGE WAVEFIELD
    figure
    hold on

    averageWavefield = mean(p.eta_data, 3);
    DrawWavefield2D(p, averageWavefield)
    
    title("Average Wavefield", 'Interpreter', 'latex')
    xlabel('$x/\lambda_F$', Interpreter = "latex")
    ylabel('$y/\lambda_F$', Interpreter = "latex")
    axis square
    xlim(BOUNDS_R_PERL_CROPPED)
    ylim(BOUNDS_R_PERL_CROPPED)
    colorbar
    clim(WAVEFIELD_CLIM)

    exportgraphics(gca, pathId + "_average_wavefield.png");

    %% WAVEFIELD CROSS X
    figure
    wavefield = p.eta_data(:, :, end);
    droplet1_x = p.x_data(end, 1);
    droplet1_y = p.y_data(end, 1);
    wavefieldX = interp2(p.xx, p.yy, wavefield, p.xx(1, :), droplet1_y);
    wavefieldDroplet = interp2(p.xx, p.yy, wavefield, droplet1_x, droplet1_y);

    plot(p.xx(1, :), wavefieldX);
    hold on
    scatter(droplet1_x, wavefieldDroplet, CROSS_SECTIONS_DROPLET_SIZE, 'filled');
    DrawBounds1DDia(p);
    hold off

    title("Wavefield (X Cross-Section)", 'Interpreter', 'latex')
    xlabel('$x/\lambda_F$', Interpreter = "latex")
    ylabel('$\eta$', Interpreter = "latex")
    xlim(BOUNDS_R_PERL)
    ylim(WAVEFIELD_CLIM)

    exportgraphics(gca, pathId + "_wavefield_x_cross.png");

    %% WAVEFIELD CROSS Y
    figure
    wavefield = p.eta_data(:, :, end);
    droplet1_x = p.x_data(end, 1);
    droplet1_y = p.y_data(end, 1);
    wavefieldY = interp2(p.xx, p.yy, wavefield, droplet1_x, p.yy(:, 1));

    plot(p.yy(:, 1), wavefieldY);
    hold on
    scatter(droplet1_y, wavefieldDroplet, CROSS_SECTIONS_DROPLET_SIZE, 'filled');
    DrawBounds1DDia(p);
    hold off

    title("Wavefield (Y Cross-Section)", 'Interpreter', 'latex')
    xlabel('$y/\lambda_F$', Interpreter = "latex")
    ylabel('$\eta$', Interpreter = "latex")
    xlim(BOUNDS_R_PERL)
    ylim(WAVEFIELD_CLIM)

    exportgraphics(gca, pathId + "_wavefield_y_cross.png");

    %% TRAJECTORIES (VIDEO)
    videoWriter = VideoWriter(pathId + "_trajectories.avi", 'Motion JPEG AVI');
    videoWriter.FrameRate = 1/p.TF;
    videoWriter.Quality = 95;
    videoWriter.open();

    figureHandle = figure(Visible = "off");

    bufferFile = tempname() + ".png";

    frameCount = size(p.x_data, 1);
    for j = 1:min(frameCount, 2 * 60 / p.TF) % 2 minutes

        % Set up shared data
        dropletPositions = zeros(p.n_drops, 2);
        for k = 1:size(p.x_data, 2)
            dropletPositions(k, :) = [p.x_data(j, k), p.y_data(j, k)];
        end

        % Wavefield
        imagesc(axisValuesX, axisValuesY, bins');
        hold on
        viscircles(dropletPositions, p.drop_radius / p.lambdaF * ones(1, p.n_drops));
        DrawBounds2D(p)
        hold off

        title("Trajectories", 'Interpreter', 'latex')
        xlabel('$x/\lambda_F$', Interpreter = "latex")
        ylabel('$y/\lambda_F$', Interpreter = "latex")

        axis square
        xlim(BOUNDS_R_PERL_CROPPED)
        ylim(BOUNDS_R_PERL_CROPPED)
        colorbar

        % Write out results
        exportgraphics(gca, bufferFile, 'Resolution', 300);

        frame = imread(bufferFile);

        videoWriter.writeVideo(frame);
    end

    videoWriter.close();

    %% WAVEFIELD (VIDEO)
    videoWriter = VideoWriter(pathId + "_wavefield.avi", 'Motion JPEG AVI');
    videoWriter.FrameRate = 1/p.TF;
    videoWriter.Quality = 95;
    videoWriter.open();

    figureHandle = figure(Visible = "off");
    layoutHandle = tiledlayout(2, 2 , Padding = "tight", TileSpacing = "tight");

    bufferFile = tempname() + ".png";

    frameCount = size(p.eta_data, 3);
    for j = 1:frameCount
           
        % Set up shared data
        absoluteTime = p.nimpacts - p.n_save_wave + j;
        dropletPositions = zeros(p.n_drops, 2);
        for k = 1:size(p.x_data, 2)
            dropletPositions(k, :) = [p.x_data(absoluteTime, k), p.y_data(absoluteTime, k)];
        end
        droplet1Position = dropletPositions(1, :);
        wavefieldX = interp2(p.xx, p.yy, wavefield, p.xx(1, :), droplet1Position(2));
        wavefieldY = interp2(p.xx, p.yy, wavefield, droplet1Position(1), p.yy(:, 1));
        wavefieldDroplet = interp2(p.xx, p.yy, wavefield, droplet1Position(1), droplet1Position(2));

        % Wavefield
        tile1Axes = nexttile(1);
        hold(tile1Axes, 'on');
        cla(tile1Axes);
        wavefield = p.eta_data(:, :, j);
        DrawWavefield2D(p, wavefield)
        viscircles(dropletPositions, p.drop_radius / p.lambdaF * ones(1, p.n_drops));
        hold(tile1Axes, 'off');

        title("Wavefield", 'Interpreter', 'latex')
        xlabel('$x/\lambda_F$', Interpreter = "latex")
        ylabel('$y/\lambda_F$', Interpreter = "latex")

        axis square
        xlim(BOUNDS_R_PERL)
        ylim(BOUNDS_R_PERL)
        colorbar
        clim(WAVEFIELD_CLIM)

        % Wavefield (Tracked)
        tile2Axes = nexttile(2);
        hold(tile2Axes, 'on');
        cla(tile2Axes);
        wavefield = p.eta_data(:, :, j);
        DrawWavefield2D(p, wavefield, droplet1Position)
        viscircles(zeros(p.n_drops, 2), p.drop_radius / p.lambdaF * ones(1, p.n_drops));
        hold(tile2Axes, 'off');

        title("Wavefield (Tracked)", 'Interpreter', 'latex')
        xlabel('$x/\lambda_F$', Interpreter = "latex")
        ylabel('$y/\lambda_F$', Interpreter = "latex")

        axis square
        xlim(BOUNDS_R_PERL)
        ylim(BOUNDS_R_PERL)
        colorbar
        clim(WAVEFIELD_CLIM)

        % Wavefield Cross X
        tile3Axes = nexttile(3);
        hold(tile3Axes, 'on');
        cla(tile3Axes);
        plot(p.xx(1, :), wavefieldX);
        scatter(droplet1Position(1), wavefieldDroplet, CROSS_SECTIONS_DROPLET_SIZE, 'filled');
        DrawBounds1DDia(p);
        hold(tile3Axes, 'off');

        title("Wavefield (X Cross-Section)", 'Interpreter', 'latex')
        xlabel('$x/\lambda_F$', Interpreter = "latex")
        ylabel('$\eta$', Interpreter = "latex")
        axis square
        xlim(BOUNDS_R_PERL)
        ylim(WAVEFIELD_CLIM)

        % Wavefield Cross Y
        tile4Axes = nexttile(4);
        hold(tile4Axes, 'on');
        cla(tile4Axes);
        plot(p.yy(:, 1), wavefieldY);
        scatter(droplet1Position(2), wavefieldDroplet, CROSS_SECTIONS_DROPLET_SIZE, 'filled');
        DrawBounds1DDia(p);
        hold(tile4Axes, 'off');

        title("Wavefield (Y Cross-Section)", 'Interpreter', 'latex')
        xlabel('$y/\lambda_F$', Interpreter = "latex")
        ylabel('$\eta$', Interpreter = "latex")
        axis square
        xlim(BOUNDS_R_PERL)
        ylim(WAVEFIELD_CLIM)

        % Write out results
        exportgraphics(layoutHandle, bufferFile, 'Resolution', 300);

        frame = imread(bufferFile);

        videoWriter.writeVideo(frame);
    end

    videoWriter.close();
end

function DrawBounds1DRad(p, origin)
    arguments
        p
        origin (1,1) double = 0
    end

    hold on

    switch p.topography_type
        case 'circular_well'
            xline(p.Rc - origin, "-", LineWidth = 1.0);
    end
    switch p.damping_type
        case 'scaled'
            xline(p.effective_corral_radius - origin, "--r", LineWidth = 1.0);
    end
    switch p.corral_type
        case 'spring'
            xline(p.effective_corral_radius - origin, "--r", LineWidth = 1.0);
        case 'rigid'
            xline(p.Rc - origin, "-", LineWidth = 1.0);
    end

    hold off
end

function DrawBounds1DDia(p, origin)
    arguments
        p
        origin (1,1) double = 0
    end

    hold on

    switch p.topography_type
        case 'circular_well'
            xline(p.Rc - origin, "-", LineWidth = 1.0);
            xline(-(p.Rc - origin), "-", LineWidth = 1.0);
    end
    switch p.damping_type
        case 'scaled'
            xline(p.effective_corral_radius - origin, "--r", LineWidth = 1.0);
            xline(-(p.effective_corral_radius - origin), "--r", LineWidth = 1.0);
    end
    switch p.corral_type
        case 'spring'
            xline(p.effective_corral_radius - origin, "--r", LineWidth = 1.0);
            xline(-(p.effective_corral_radius - origin), "--r", LineWidth = 1.0);
        case 'rigid'
            xline(p.Rc - origin, "-", LineWidth = 1.0);
            xline(-(p.Rc - origin), "-", LineWidth = 1.0);
    end
    
    hold off
end

function DrawWavefield2D(p, wavefield, origin)
    arguments
        p
        wavefield
        origin (1,2) double = [0 0]
    end

    contourf(p.xx  + p.Lx/p.Nx/2 - origin(1), p.yy + p.Ly/p.Ny/2 - origin(2), wavefield, 50, "EdgeColor", "none");
    hold on
    DrawBounds2D(p, origin)
    hold off
end

function DrawBounds2D(p, origin)
    arguments
        p
        origin (1,2) double = [0 0]
    end

    hold on
    
    viscircles([0 0] - origin, [p.Rc], LineWidth = 0.2,  LineStyle = '-', Color = 'black');
    switch p.damping_type
        case 'scaled'
            viscircles([0 0] - origin, [p.effective_corral_radius], LineWidth = 0.1,  LineStyle = '--', Color = 'red');
    end
    switch p.corral_type
        case 'spring'
            viscircles([0 0] - origin, [p.effective_corral_radius], LineWidth = 0.1,  LineStyle = '--', Color = 'red');
    end
    
    hold off
end