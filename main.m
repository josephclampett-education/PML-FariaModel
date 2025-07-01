%% Set Fluid Parameters

p = fluid_params();

%% Set Topography Parameters

% Domain Size
p.Lx = 19; % lambdaF
p.Ly = 19; % lambdaF

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
        p.h0 = 1.5*10^(-3); % m (constant depth)
    case 'square_well'
        p.h0 = 1.5*10^(-3); % m (interior depth)
        p.h1 = 1.5*10^(-4); % m (exterior depth)
        p.Lt = 4;           % lambdaF (well width)
    case 'circular_well'
        p.h0 = 1.5*10^(-3); % m (interior depth)
        p.h1 = 1.5*10^(-4); % m (exterior depth)
        p.Dt = 4;           % lambdaF (well diameter)
        p.Rt = p.Dt/2;      % lambdaF (well radius)
end

p = top_params(p);

%% Calculate Faraday Threshold

p = faraday_threshold(p);

%% Set Drop Parameters & Initial Conditions

% Drop Size
p.drop_radius  = (0.745/2)*10^(-3);                       % m
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
p.n_drops = 1;

% Set Drop Initial Conditions
p.xi = zeros(1,p.n_drops);
p.yi = zeros(1,p.n_drops);
p.ui = zeros(1,p.n_drops)+0.01;
p.vi = zeros(1,p.n_drops);

% Set Wave Initial Conditions
p.eta0 = zeros(size(p.xx));
p.phi0 = zeros(size(p.xx));

% Set Memory
p.mem = 0.95;
p.Gam = p.mem*p.GamF;

% Set Number of Impacts (Simulation Time in TF)
p.nimpacts = 500;

p = drop_params(p);

%% Run Simulation

p = simulate(p);