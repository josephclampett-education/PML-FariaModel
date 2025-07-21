% ================================================================
% REMINDERS
% ================================================================

% Make sure to always check that the core loop is parfor not for!
% If making a new run id (ex: 09 or 51) make sure to notify one on Slack

% ================================================================
% PARAMS
% ================================================================

BASE_DIRECTORY = "../..";

addpath(BASE_DIRECTORY);

% BATCH
BATCH0_h1 = [0.2, 0.3, 0.5, 0.8] * 10^(-3);
BATCH0_h0 = BATCH0_h1 + (5.46*10^(-3) - 0.61*10^(-3));
BATCH1_R = linspace(2.37633, 2.8761, 5);

% BATH
VAR_mem = 0.99;
VAR_type = 'circular_well';
% VAR_R = 4.8757;
% VAR_h0 = 5.46*10^(-3);
% VAR_h1 = 0.61*10^(-3);

VAR_shouldOverrideThreshold = 0;
VAR_thresholdGuess = 4.1979;

% DROPLETS
VAR_r = (0.36)*10^(-3);
VAR_sin_theta = 0.2;
VAR_n_drops = 10;

% INITIAL CONDITIONS
VAR_initialRadiusScale = 0.80;
VAR_initialSpeedScale = 0.01;

% SIMULATION
if isfile(BASE_DIRECTORY + "/ISLOCAL")
    VAR_nimpacts = 2;
    VAR_n_save_wave = 1;
else
    VAR_nimpacts = 2;
    VAR_n_save_wave = 1;
end

% Saving
VAR_outputFolder = "RES";

%% ================================================================

count0 = length(BATCH0_h1);
count1 = length(BATCH1_R);
threadCount = count0 * count1;

outputData = [];

parfor i = 1:threadCount
    %% Unpack Dispatch Parameters
    idx0 = mod((i - 1), count0) + 1;
    idx1 = floor((i - 1) / count0) + 1;

    %% Set Fluid Parameters
    
    p = fluid_params();
    
    %% Set Topography Parameters
    
    % Domain Size
    p.Lx = 2 * 8; % lambdaF
    p.Ly = 2 * 8; % lambdaF
    
    % Grid Spacing
    p.hx_desired = 0.1;                          % lambdaF
    p.hy_desired = 0.1;                          % lambdaF
    p.hx         = p.Lx/ceil(p.Lx/p.hx_desired); % lambdaF
    p.hy         = p.Ly/ceil(p.Ly/p.hy_desired); % lambdaF
    
    p.Nx = p.Lx/p.hx; % dimensionless
    p.Ny = p.Ly/p.hy; % dimensionless
    
      % Note:
      % round desired grid spacing to whole number of points per lambdaF
    
    % Time Step (for wave evolution)
    p.dt_desired = min(p.hx,p.hy)^2/10;    % TF
    p.dt         = 1/ceil(1/p.dt_desired); % TF
    
    p.nsteps_impact = 1/p.dt; % dimensionless
    
      % Note:
      % currently set based on spatial resolution, could set independently
      % does not affect drop since instantaneous impacts, only wave evolution
      % round desired time step to have whole number of steps per TF
      % scale with spatial res squared to keep well behaved for high res
    
    % Topography
    p.type = VAR_type; % options: 'flat', 'square_well', 'circular_well'
    
    radius = BATCH1_R(idx1);
    switch p.type
        case 'flat'
            p.h0 = BATCH0_h0(idx0); % m (constant depth)
            p.h1 = p.h0;
            p.Rc = radius;  % lambdaF (droplet corral radius)
            p.Dc = p.Rc*2; % lambdaF (droplet corral diameter)
        case 'square_well'
            p.h0 = 1.5*10^(-3); % m (interior depth)
            p.h1 = 1.5*10^(-4); % m (exterior depth)
            p.Lt = 4;           % lambdaF (well width)
        case 'circular_well'
            p.h0 = BATCH0_h0(idx0); % m (interior depth)
            p.h1 = BATCH0_h1(idx0); % m (exterior depth)
            p.Rc = radius;  % lambdaF (well radius)
            p.Dc = p.Rc*2; % lambdaF (well diameter)
    end
    
    p = top_params(p);
    
    %% Calculate Faraday Threshold
    if VAR_shouldOverrideThreshold
        fprintf("%s: Overriding threshold.\n", datetime);
        p.GamF = VAR_thresholdGuess;
    else
        thresholdFile = sprintf("%s/threshold_cache/%f_%f_%f.mat", BASE_DIRECTORY, p.h0, p.h1, p.Rc);
        if isfile(thresholdFile)
          fprintf("%s: Cache hit, loading threshold for %s.\n", datetime, thresholdFile);
          gamFLoad = load(thresholdFile);
          p.GamF = gamFLoad.GamF;
        else
          fprintf("%s: Cache miss, calculating threshold for %s.\n", datetime, thresholdFile);
          p = faraday_threshold(p, VAR_thresholdGuess);
          GamF = p.GamF;
          parsave_gamf(thresholdFile, GamF);
        end
    end
    
    %% Set Drop Parameters & Initial Conditions
    
    % Drop Size
    p.drop_radius  = VAR_r;                                   % m
    p.drop_density = 949;                                     % kg/m^3
    p.drop_mass    = (4/3)*pi*p.drop_radius^3*p.drop_density; % kg
    
      % Note:
      % only the mass matters since treated as a point for impacts
    
    % Impact Phase
    p.sin_theta = VAR_sin_theta;
    p.theta     = pi-asin(p.sin_theta);
    
      % Note:
      % effectively controls speed of drop given other parameters
      % optional: can choose a phase to match speed shown in experiments
    
    % Number of Drops
    p.n_drops = VAR_n_drops;
    
    % Set Drop Initial Conditions
    randTheta = 2*pi*rand(1,p.n_drops);
    adjustedMaxRandR = VAR_initialRadiusScale * p.Rc;
    randR = adjustedMaxRandR*sqrt(rand(1,p.n_drops));

    randX = cos(randTheta);
    randY = sin(randTheta);

    p.xi = randR .* randX;
    p.yi = randR .* randY;
    p.ui = VAR_initialSpeedScale * rand(1,p.n_drops);
    p.vi = VAR_initialSpeedScale * rand(1,p.n_drops);
    
    % Set Wave Initial Conditions
    p.eta0 = zeros(size(p.xx));
    p.phi0 = zeros(size(p.xx));
    
    % Set Memory
    p.mem = VAR_mem;
    p.Gam = p.mem*p.GamF;
    
    % Set Number of Impacts (Simulation Time in TF)
    p.nimpacts = VAR_nimpacts;
    
    p = drop_params(p);

    % Save Parameters
    p.n_save_wave = VAR_n_save_wave;
    
    %% Run Simulation
    
    fprintf("%s: Beginning simulation.\n", datetime);
    p = simulate(p);

    %% Output Results
    if ~isfolder(VAR_outputFolder)
        mkdir(VAR_outputFolder);
    end
    saveFilePath = sprintf("%s/RES_mem=%f, N=%d, %s R=%f.mat", VAR_outputFolder, p.mem, p.n_drops, p.type, p.Rc);
    fprintf("%s: Saving simulation results for %s.\n", datetime, saveFilePath);
    parsave(saveFilePath, p);
end

%% Hack to allow saving inside parfor
function parsave(fname, p)
  save(fname, 'p', "-v7.3");
end

%% Hack to allow saving inside parfor
function parsave_gamf(fname, GamF)
  save(fname, 'GamF', "-v7.3");
end

