addpath("..");

memoryVals = [0.95, 0.98, 0.99, 0.995];
countVals = 10:10:40;

count0 = length(memoryVals);
count1 = length(countVals);
threadCount = count0 * count1;

outputData = [];

parfor i = 1:threadCount
    %% Unpack Dispatch Parameters
    idx0 = mod((i - 1), count0) + 1;
    idx1 = floor((i - 1) / count0) + 1;

    mem = memoryVals(idx0);
    count = countVals(idx1);

    %% Set Fluid Parameters
    
    p = fluid_params();
    
    %% Set Topography Parameters
    
    % Domain Size
    p.Lx = 20; % lambdaF
    p.Ly = 20; % lambdaF
    
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
    p.type = 'flat'; % options: 'flat', 'square_well', 'circular_well'
    
    switch p.type
        case 'flat'
            p.h0 = 4*10^(-3); % m (constant depth)
            p.Dp = 4.8757*2;  % lambdaF (droplet corral diameter)
            p.Rp = p.Dp/2;    % lambdaF (droplet corral radius)
        case 'square_well'
            p.h0 = 1.5*10^(-3); % m (interior depth)
            p.h1 = 1.5*10^(-4); % m (exterior depth)
            p.Lt = 4;           % lambdaF (well width)
        case 'circular_well'
            % JX -- In the flat case we're using Ian's version, in this case
            % we're using experimental parameters for depths but using Ian's
            % case for circular radius
            p.h0 = 5.46*10^(-3); % m (interior depth)
            p.h1 = 0.61*10^(-3); % m (exterior depth)
            p.Dt = 4.8757*2;     % lambdaF (well diameter)
            p.Rt = p.Dt/2;       % lambdaF (well radius)
    end
    
    p = top_params(p);
    
    %% Calculate Faraday Threshold
    
    % p = faraday_threshold(p);
    p.nF = 1000;
    p.GamF = 4.1979;
    
    %% Set Drop Parameters & Initial Conditions
    
    % Drop Size
    % JX -- TODO, Need to reconcile this size with the strobe model
    p.drop_radius  = (0.36)*10^(-3);                          % m
    p.drop_density = 949;                                     % kg/m^3
    p.drop_mass    = (4/3)*pi*p.drop_radius^3*p.drop_density; % kg
    
      % Note:
      % only the mass matters since treated as a point for impacts
    
    % Impact Phase
    p.sin_theta = 0.2;
    p.theta     = pi-asin(p.sin_theta);
    
      % Note:
      % effectively controls speed of drop given other parameters
      % optional: can choose a phase to match speed shown in experiments
    
    % Number of Drops
    p.n_drops = count;
    
    % Set Drop Initial Conditions
    randTheta = 2*pi*rand(1,p.n_drops);
    adjustedMaxRandR = 0.80 * p.Rp;
    randR = adjustedMaxRandR*sqrt(rand(1,p.n_drops));

    randX = cos(randTheta);
    randY = sin(randTheta);

    p.xi = randR .* randX;
    p.yi = randR .* randY;
    p.ui = 0.01 * randX;
    p.vi = 0.01 * randY;
    
    % Set Wave Initial Conditions
    p.eta0 = zeros(size(p.xx));
    p.phi0 = zeros(size(p.xx));
    
    % Set Memory
    p.mem = mem;
    p.Gam = p.mem*p.GamF;
    
    % Set Number of Impacts (Simulation Time in TF)
    p.nimpacts = 40 * 60 * 5;
    
    p = drop_params(p);
    
    %% Run Simulation
    
    p = simulate(p);

    %% Output Results
    saveFilePath = fullfile(runFolder, ['mem=', num2str(p.mem), ' ', 'N=', num2str(p.n_drops), '.mat']);
    parsave(saveFilePath, p);
end

%% Hack to allow saving inside parfor
function parsave(fname, p)
  save(fname, 'p');
end
