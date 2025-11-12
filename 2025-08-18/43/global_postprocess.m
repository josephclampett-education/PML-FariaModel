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

POST_memories = [];
POST_phases = [];
POST_memoryVelocityPairs = [];

POST_memoryEtaPairs = [];

parfor i = 1:threadCount
    folderPath = fullfile(folders(i).folder, folders(i).name);

    % pathId is important for saving!
    [pathstr, name, ext] = fileparts(folderPath);
    pathId = "GLOBALRESFIG";

    fileSearchPath = fullfile(folderPath, "*.mat");
    files = dir(fileSearchPath);

    filePath = fullfile(files(1).folder, files(1).name);
    loadP = load(filePath);
    p = loadP.p;

    close all

    %% Setup

    bounds = [-(p.Rc + 1), (p.Rc + 1)];

    memory = p.mem;
    phase = p.theta / pi;

    if (~any(POST_memories(:) == memory))
        POST_memories = [POST_memories, memory];
    end
    if (~any(POST_phases(:) == phase))
        POST_phases = [POST_phases, phase];
    end
    
    %% VELOCITY POSTPROCESS

    xs_pl = p.x_data;
    ys_pl = p.y_data;

    vx_plpf = diff(xs_pl, 1, 1);
    vx_plpf = mod(vx_plpf + p.Lx/2, p.Lx) - p.Lx/2;

    vy_plpf = diff(ys_pl, 1, 1);
    vy_plpf = mod(vy_plpf + p.Ly/2, p.Ly) - p.Ly/2;

    vs_mmpf = sqrt(vx_plpf.^2 + vy_plpf.^2) * p.lambdaF * 1000;
    vs_mmps = vs_mmpf * 1/p.TF;

    % Get average velocity of last batch
    averageVelocity = mean(vs_mmps(floor(end/2):end));
    POST_memoryVelocityPairs = [POST_memoryVelocityPairs; [memory, averageVelocity]];

    %% WAVEFIELD POSTPROCESS

    wavefield = p.eta_data(:, :, end);
    droplet1_x = p.x_data(end, 1);
    droplet1_y = p.y_data(end, 1);
    droplet1_eta = interp2(p.xx, p.yy, wavefield, droplet1_x, droplet1_y);

    POST_memoryEtaPairs = [POST_memoryEtaPairs; [memory, droplet1_eta]];
end

%% VELOCITY VS PHASE
figure

hold on
for i = 1:size(POST_memories, 2)
    memory = POST_memories(i);

    filteredVelocities = [];
    for j = 1:size(POST_memoryVelocityPairs, 1)
        pair = POST_memoryVelocityPairs(j, :);
        if pair(1) == memory
            filteredVelocities = [filteredVelocities, pair(2)];
        end
    end

    scatter(POST_phases, filteredVelocities, 'o', MarkerEdgeColor = "black", HandleVisibility = "off");
    plot(POST_phases, filteredVelocities, DisplayName = sprintf("%.0f%%", memory*100));
end
hold off

legend;

title("Velocity vs. Phase", 'Interpreter', 'latex')
xlabel('$\phi$', Interpreter = "latex")
ylabel('$v$ (mm/s)', Interpreter = "latex")
oldYlim = ylim;
ylim([-1.0, oldYlim(2)])

exportgraphics(gca, fullfile(VAR_outputFolder, pathId + "_velocity_vs_phase.png"));

%% WAVEFIELD VS PHASE
figure

hold on
for i = 1:size(POST_memories, 2)
    memory = POST_memories(i);

    filteredEtas = [];
    for j = 1:size(POST_memoryEtaPairs, 1)
        pair = POST_memoryEtaPairs(j, :);
        if pair(1) == memory
            filteredEtas = [filteredEtas, pair(2)];
        end
    end

    scatter(POST_phases, filteredEtas, 'o', MarkerEdgeColor = "black", HandleVisibility = "off");
    plot(POST_phases, filteredEtas, DisplayName = sprintf("%.0f%%", memory*100));
end
hold off

lgd = legend;

title("Wavefield vs. Phase", Interpreter = "latex")
xlabel('$\phi$', Interpreter = "latex")
ylabel('$\eta$ (m)', Interpreter = "latex")

exportgraphics(gca, fullfile(VAR_outputFolder, pathId + "_wavefield_vs_phase.png"));