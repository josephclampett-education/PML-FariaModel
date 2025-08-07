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

% <<< MATT >>> I've changed the value of the impact phase to account for
% the change in definition of vibrational acceleration below. The new value
% of theta might need to be tweaked, but it should be about the value that
% I've listed below. I expect that we will have 1.35 <= theta/pi <= 1.55.
% The closer theta gets to pi, the faster the droplet walks. The further
% theta gets from pi, the slower the droplet walks. I used theta = 1.43 *
% pi in the corral paper, for example.

% BATCH
BATCH_h1 = [0.1] * 10^(-3);
BATCH_theta = [1.20];

% BATH
VAR_type = 'circular_well';

VAR_h0_base = 4.85*10^(-3);    % mm
% VAR_h1 = 0.20*10^(-3);    % mm
VAR_R  = 2.376330;          % in xF

VAR_mem = 0.95;

VAR_shouldOverrideThreshold = 0;
VAR_thresholdGuess = 5.0166;

% DROPLETS
VAR_r = (0.36)*10^(-3);
% VAR_theta = 1.3;
VAR_n_drops = 1;

% INITIAL CONDITIONS
VAR_initialRadiusScale = 0.80;
VAR_initialSpeedScale = 0.01;

% SIMULATION
if isfile(BASE_DIRECTORY + "/ISLOCAL")
    VAR_nimpacts = 10;
    VAR_n_save_wave = 10;
else
    VAR_nimpacts = 40*60;
    VAR_n_save_wave = 10;
end

% Saving
VAR_outputFolder = "RES";

%% ================================================================

count0 = length(BATCH_h1);
count1 = length(BATCH_theta);
threadCount = count0 * count1;

% Only do one run if using on local
if isfile(BASE_DIRECTORY + "/ISLOCAL")
  threadCount = 1;
end

outputData = [];

for i = 1:threadCount
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
    p.hx_desired = p.Lx/256;                          % lambdaF
    p.hy_desired = p.Ly/256;                          % lambdaF
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
    
    radius = VAR_R;
    switch p.type
        case 'flat'
            p.h1 = VAR_h0_base;
            p.h0 = VAR_h0_base; % m (constant depth)
            p.Rc = radius; % lambdaF (droplet corral radius)
            p.Dc = p.Rc*2; % lambdaF (droplet corral diameter)
        case 'square_well'
            p.h0 = 1.5*10^(-3); % m (interior depth)
            p.h1 = 1.5*10^(-4); % m (exterior depth)
            p.Lt = 4;           % lambdaF (well width)
        case 'circular_well'
            p.h1 = BATCH_h1(idx0);     % m (interior depth)
            p.h0 = VAR_h0_base + p.h1; % m (exterior depth)
            p.Rc = radius;  % lambdaF (well radius)
            p.Dc = p.Rc*2;  % lambdaF (well diameter)
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
end

%% Hack to allow saving inside parfor
function parsave(fname, p)
  save(fname, 'p', "-v7.3");
end

%% Hack to allow saving inside parfor
function parsave_gamf(fname, GamF)
  save(fname, 'GamF', "-v7.3");
end

