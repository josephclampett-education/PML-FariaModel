addpath("..");

memoryVals = [0.95, 0.98, 0.99, 0.995];
countVals = 10:10:20;

counts = [length(memoryVals), length(countVals)];
threadCount = prod(counts);

parfor i = 1:threadCount
    %% Unpack Dispatch Parameters
    [idx0, idx1] = ind2sub(counts, i);

    mem = memoryVals(idx0);
    count = countVals(idx1);

    bathData = load("bath.mat");
    p = bathData.p;

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
    adjustedMaxRandR = 0.80 * p.R0;
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

    % Save Parameters
    p.n_save_wave = 10; % Save wavefield for the last n impacts
    
    %% Run Simulation
    p = simulate(p);

    %% Output Results
    saveFilePath = fullfile(runFolder, ['mem=', num2str(p.mem), ' ', 'N=', num2str(p.n_drops), '.mat']);
    parsave(saveFilePath, p);
end

%% Hack to allow saving inside parfor
function parsave(fname, p)
  save(fname, 'p', "-v7.3");
end
