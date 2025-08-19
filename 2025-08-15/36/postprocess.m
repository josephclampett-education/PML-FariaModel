BASE_DIRECTORY = "../..";

addpath(BASE_DIRECTORY);

% ================================================================
% PARAMS
% ================================================================

% Saving
VAR_outputFolder = "RES";

%% ================================================================

folders = dir(fullfile(VAR_outputFolder, "RES_*"));
threadCount = length(folders);

outputData = [];

threadCount = 1;
for i = 1:threadCount
    folderPath = fullfile(folders(i).folder, folders(i).name);

    % pathId is important for saving!
    [pathstr, name, ext] = fileparts(folderPath);
    pathId = fullfile(folderPath, strcat(name, ext)); % Add duplicate of innermost folder string

    fileSearchPath = fullfile(folderPath, "*.mat");
    files = dir(fileSearchPath);

    x_data = [];
    y_data = [];
    for j = 1:numel(files)
        filePath = fullfile(files(j).folder, files(j).name);
        loadP = load(filePath);
        p = loadP.p;

        x_data = [x_data; p.x_data];
        y_data = [y_data; p.y_data];
    end

    close all

    %% Setup
    xs = [];
    ys = [];
    for j = 1:size(x_data, 2)
        xs = [xs; x_data(:, j)];
        ys = [ys; y_data(:, j)];
    end

    bounds = [-(p.Rc + 1), (p.Rc + 1)];
    
    %% HISTOGRAM
    figure
    cellcount = 180;
    
    xEdges = linspace(-p.Lx/2, p.Lx/2, cellcount + 1);
    yEdges = linspace(-p.Ly/2, p.Ly/2, cellcount + 1);
    [bins, xEdges, yEdges] = histcounts2(xs, ys, xEdges, yEdges, Normalization="probability");

    axisValuesX = linspace(-p.Lx/2, p.Lx/2, cellcount);
    axisValuesY = linspace(-p.Ly/2, p.Ly/2, cellcount);
    hold on
    imagesc(axisValuesX, axisValuesY, bins');
    viscircles([0, 0], [p.Rc], 'LineWidth', 0.2, 'LineStyle','--');
    hold off
    
    axis square
    %colorbar
    xlim(bounds)
    ylim(bounds)
    %clim([0.1, 1.7] * 10e-5)

    exportgraphics(gca, pathId + "_histogram.png");
    
    %% RADIAL HISTOGRAM
    figure
    bincount = 80;
    
    r_data = sqrt(xs.^2 + ys.^2);
    radEdges = linspace(0, p.Rc, bincount);
    radBins = histcounts(r_data, radEdges);
    radCenters = (radEdges(1:end - 1) + radEdges(2:end)) / 2;
    radProb = radBins ./ radCenters;
    radProb = radProb / length(r_data);
    
    j0Domain = (0:0.02:p.Rc);
    J0 = besselj(0, (2*pi) * j0Domain);
    
    hold on
    bar(radCenters, radProb, 'hist');
    plot(j0Domain, abs(J0) * max(radProb), Color="red");
    xlim([0 p.Rc]);
    hold off

    exportgraphics(gca, pathId + "_radial_histogram.png");

    %% TRAJECTORY
    figure
    hold on
    for j = 1:size(p.x_data, 2)
        plot(x_data(:, j), y_data(:, j));
    end
    viscircles([0, 0], [p.Rc], 'LineWidth', 0.2, 'LineStyle','--');
    hold off
    
    axis square
    xlim(bounds)
    ylim(bounds)

    exportgraphics(gca, pathId + "_trajectory.png");

    %% WAVEFIELD
    figure
    hold on
    wavefield = p.eta_data(:, :, end);
    contourf(p.xx, p.yy, wavefield, 50, "EdgeColor", "none");
    viscircles([p.x_data(end,:); p.y_data(end,:)]', p.drop_radius / p.lambdaF * ones(1, p.n_drops));
    viscircles([0, 0], [p.Rc], 'LineWidth', 0.2, 'LineStyle','--');
    hold off
    
    axis square
    xlim(bounds)
    ylim(bounds)

    exportgraphics(gca, pathId + "_wavefield.png");

    %% WAVEFIELD (VIDEO)
    v = VideoWriter(pathId + "_wavefield.mp4", 'MPEG-4');
    v.FrameRate = 1/p.TF;
    open(v);

    fig5 = figure(5);
    frameCount = size(p.eta_data, 3);
    for j = 1:frameCount
        % Wavefield
        wavefield = p.eta_data(:, :, j);
        contourf(p.xx, p.yy, wavefield, 50, "EdgeColor", "none");

        hold on
        % Droplet positions
        absoluteTime = p.nimpacts - p.n_save_wave + j;
        dropletPositions = zeros(p.n_drops, 2);
        for k = 1:size(p.x_data, 2)
            dropletPositions = [p.x_data(absoluteTime, k), p.y_data(absoluteTime, k)];
        end
        viscircles(dropletPositions, p.drop_radius / p.lambdaF * ones(1, p.n_drops));

        % Corral
        viscircles([0, 0], [p.Rc], 'LineWidth', 0.2, 'LineStyle','--');
        hold off

        axis square
        xlim(bounds)
        ylim(bounds)
        title(sprintf("Wavefield at timestep %d / %d", j, frameCount));

        frame = getframe(gcf);
        writeVideo(v, frame);
    end

    close(v);
end