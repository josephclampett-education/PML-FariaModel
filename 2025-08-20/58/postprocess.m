BASE_DIRECTORY = "../..";

addpath(BASE_DIRECTORY);

% ================================================================
% PARAMS
% ================================================================

HISTOGRAM_BINCOUNT = 180;
HISTOGRAM_CLIM = [0.1 1.1] * 10^-3;

RADIAL_HISTOGRAM_BINCOUNT = 80;
RADIAL_HISTOGRAM_USEPLIM = false;
RADIAL_HISTOGRAM_PLIM = 0.025;

VELOCITY_HISTOGRAM_BINCOUNT = RADIAL_HISTOGRAM_BINCOUNT;
VELOCITY_HISTOGRAM_MAXVEL = 40;

WAVEFIELD_CLIM = [-1 +1] * 0.010;

CROSS_SECTIONS_DROPLET_SIZE = 20;

% Saving
VAR_outputFolder = "RES";

%% ================================================================

folders = dir(fullfile(VAR_outputFolder, "RES_*"));
threadCount = length(folders);

% % Only do one run if using on local
% if isfile(BASE_DIRECTORY + "/ISLOCAL")
%   threadCount = 1;
% end

for i = 1:threadCount
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
    hold on
    for j = 1:size(p.x_data, 2)
        plot(x_data(:, j), y_data(:, j));
    end
    DrawBounds2D(p)

    title("Trajectories", 'Interpreter', 'latex')
    xlabel('$x/\lambda_F$', Interpreter = "latex")
    ylabel('$y/\lambda_F$', Interpreter = "latex")
    axis square
    xlim(BOUNDS_R_PERL_CROPPED)
    ylim(BOUNDS_R_PERL_CROPPED)

    exportgraphics(gca, pathId + "_trajectory.png");
    
    %% VELOCITY HISTOGRAM
    figure

    velEdges = linspace(0, VELOCITY_HISTOGRAM_MAXVEL, VELOCITY_HISTOGRAM_BINCOUNT + 1);
    velBins = zeros(1, VELOCITY_HISTOGRAM_BINCOUNT);
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

    title("Velocity Histogram", 'Interpreter', 'latex')
    xlabel('$v$ (mm/s)', Interpreter = "latex")
    ylabel('$p$', Interpreter = "latex")
    xlim([0 VELOCITY_HISTOGRAM_MAXVEL])

    exportgraphics(gca, pathId + "_velocity_histogram.png");
    
    %% VELOCITY PLOT
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

    title("Velocity", 'Interpreter', 'latex')
    xlabel('$t_n$', Interpreter = "latex")
    ylabel('$v$ (mm/s)', Interpreter = "latex")

    exportgraphics(gca, pathId + "_velocity_plot.png");

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
    ylabel('$\eta$ (m)', Interpreter = "latex")
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
    ylabel('$\eta$ (m)', Interpreter = "latex")
    xlim(BOUNDS_R_PERL)
    ylim(WAVEFIELD_CLIM)

    exportgraphics(gca, pathId + "_wavefield_y_cross.png");

    %% WAVEFIELD (VIDEO)
    v = VideoWriter(pathId + "_wavefield.avi", 'Motion JPEG AVI');
    v.FrameRate = 1/p.TF;
    v.Quality = 95;
    open(v);

    figureHandle = figure;
    tiledlayout(2, 2 , Padding = "tight", TileSpacing = "tight")

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
        nexttile(1)
        wavefield = p.eta_data(:, :, j);
        DrawWavefield2D(p, wavefield)
        hold on
        viscircles(dropletPositions, p.drop_radius / p.lambdaF * ones(1, p.n_drops));
        hold off

        title("Wavefield", 'Interpreter', 'latex')
        xlabel('$x/\lambda_F$', Interpreter = "latex")
        ylabel('$y/\lambda_F$', Interpreter = "latex")

        axis square
        xlim(BOUNDS_R_PERL)
        ylim(BOUNDS_R_PERL)
        colorbar
        clim(WAVEFIELD_CLIM)

        % Wavefield (Tracked)
        nexttile(2)
        wavefield = p.eta_data(:, :, j);
        DrawWavefield2D(p, wavefield, droplet1Position)
        hold on
        viscircles([0 0], p.drop_radius / p.lambdaF * ones(1, p.n_drops));
        hold off

        title("Wavefield (Tracked)", 'Interpreter', 'latex')
        xlabel('$x/\lambda_F$', Interpreter = "latex")
        ylabel('$y/\lambda_F$', Interpreter = "latex")

        axis square
        xlim(BOUNDS_R_PERL)
        ylim(BOUNDS_R_PERL)
        colorbar
        clim(WAVEFIELD_CLIM)

        % Wavefield Cross X
        nexttile(3)
        plot(p.xx(1, :), wavefieldX);
        hold on
        scatter(droplet1Position(1), wavefieldDroplet, CROSS_SECTIONS_DROPLET_SIZE, 'filled');
        DrawBounds1DDia(p);
        hold off

        title("Wavefield (X Cross-Section)", 'Interpreter', 'latex')
        xlabel('$x/\lambda_F$', Interpreter = "latex")
        ylabel('$\eta$ (m)', Interpreter = "latex")
        axis square
        xlim(BOUNDS_R_PERL)
        ylim(WAVEFIELD_CLIM)

        % Wavefield Cross Y
        nexttile(4)
        plot(p.yy(:, 1), wavefieldY);
        hold on
        scatter(droplet1Position(2), wavefieldDroplet, CROSS_SECTIONS_DROPLET_SIZE, 'filled');
        DrawBounds1DDia(p);
        hold off

        title("Wavefield (Y Cross-Section)", 'Interpreter', 'latex')
        xlabel('$y/\lambda_F$', Interpreter = "latex")
        ylabel('$\eta$ (m)', Interpreter = "latex")
        axis square
        xlim(BOUNDS_R_PERL)
        ylim(WAVEFIELD_CLIM)

        % Write out results
        frame = getframe(figureHandle);
        writeVideo(v, frame);
    end

    close(v);
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