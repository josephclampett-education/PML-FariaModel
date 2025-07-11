addpath("..");

h1 = 0.9;
corralRadius = 4.8757;

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
p.type = 'circular_well'; % options: 'flat', 'square_well', 'circular_well'
p.R0 = corralRadius;

switch p.type
    case 'flat'
        p.h0 = 4*10^(-3); % m (constant depth)
        p.Dp = p.R0*2;  % lambdaF (droplet corral diameter)
        p.Rp = p.Dp/2;    % lambdaF (droplet corral radius)
    case 'square_well'
        p.h0 = 1.5*10^(-3); % m (interior depth)
        p.h1 = 1.5*10^(-4); % m (exterior depth)
        p.Lt = 4;           % lambdaF (well width)
    case 'circular_well'
        % JX -- In the flat case we're using Ian's version, in this case
        % we're using experimental parameters for depths but using Ian's
        % case for circular radius
        p.h0 = 4*10^(-3); % m (interior depth)
        p.h1 = h1*10^(-3); % m (exterior depth)
        p.Dt = p.R0*2;     % lambdaF (well diameter)
        p.Rt = p.Dt/2;       % lambdaF (well radius)
end

p = top_params(p);

%% Calculate Faraday Threshold

p = faraday_threshold(p);

%% Output Results
save("bath.mat", "p");