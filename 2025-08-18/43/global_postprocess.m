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

for i = 1:threadCount
    folderPath = fullfile(folders(i).folder, folders(i).name);

    % pathId is important for saving!
    [pathstr, name, ext] = fileparts(folderPath);
    pathId = fullfile(folderPath, strcat(name, ext)); % Add duplicate of innermost folder string

    fileSearchPath = fullfile(folderPath, "*.mat");
    files = dir(fileSearchPath);

    filePath = fullfile(files(1).folder, files(1).name);
    loadP = load(filePath);
    p = loadP.p;

    close all

    %% Setup

    bounds = [-(p.Rc + 1), (p.Rc + 1)];
    
    %% VELOCITY POSTPROCESS

    xs_pl = p.x_data;
    ys_pl = p.y_data;

    vx_plpf = diff(xs_pl, 1, 1);
    vy_plpf = diff(ys_pl, 1, 1);

    vs_mmpf = sqrt(vx_plpf.^2 + vy_plpf.^2) * p.lambdaF * 1000;
    vs_mmps = vs_mmpf * 1/p.TF;

    % Get average velocity of last batch
    memory = p.mem;
    phase = p.theta / pi;
    averageVelocity = mean(vs_mmps(floor(end/2):end));
    POST_memoryVelocityPairs = [POST_memoryVelocityPairs; [memory, averageVelocity]];
    if (~any(POST_memories(:) == memory))
        POST_memories = [POST_memories, memory];
    end
    if (~any(POST_phases(:) == phase))
        POST_phases = [POST_phases, phase];
    end
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

    plot(POST_phases, filteredVelocities, DisplayName = sprintf("%.0f%%", memory*100));
end
hold off

lgd = legend;

title("Velocity vs. Phase", 'Interpreter', 'latex')
xlabel('$\phi$','Interpreter','latex')
ylabel('$v (mm/s)$','Interpreter','latex')
oldYlim = ylim;
ylim([-1.0, oldYlim(2)])

exportgraphics(gca, fullfile(VAR_outputFolder, "RES_GLOBAL_velocity_vs_phase.png"));