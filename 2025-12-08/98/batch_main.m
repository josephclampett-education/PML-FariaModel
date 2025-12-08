function batch_main(IN_batchIndex)
	BASE_DIRECTORY = "../..";
	addpath(BASE_DIRECTORY);

	% ================================================================
	% REMINDERS
	% ================================================================

	% Make sure to always check that the core loop is parfor not for!
	% If making a new run id (ex: 09 or 51) make sure to notify one on Slack

	% ================================================================
	% CONSTS
	% ================================================================

	CONST_RLIST = [2.37633, (2.37633 + 2.62138) / 2, 2.62138, (2.62138 + 2.87610) / 2, 2.87610];

	% ================================================================
	% PARAMS
	% ================================================================

	% BATH
	VAR_topography_type = "circular_well";
	VAR_damping_type = "none";
	VAR_corral_type = "none";

	VAR_h0_base = 4.85*10^(-3); % mm
	VAR_h1 = 0.30*10^(-3);      % mm
	VAR_R  = CONST_RLIST(1);    % in xF

	VAR_mem = 0.99;

    VAR_wave_damping_scale = [1, 5, 20];
    VAR_wave_damping_blend_rate = 0.2;

	VAR_shouldOverrideThreshold = false;
	VAR_thresholdGuess = 5;

	% DROPLETS
	VAR_r = (0.45)*10^(-3); % m
	VAR_theta = 1.20;       % / π
	VAR_n_drops = 10;

	VAR_droplet_collision_type = "none";

	VAR_c4 = 0.17;

	% INITIAL CONDITIONS
	VAR_initialRadiusScale = 0.80;
	VAR_initialSpeedScale = 0.1;

	% SIMULATION
	VAR_domainWidth = 2 * 8;
	VAR_gridPerWave = 8;
	VAR_timeStepDivisor = 10;
	if isfile(BASE_DIRECTORY + "/ISLOCAL")
		VAR_nimpacts = 200;
		VAR_n_save_wave = 200;
	else
		VAR_nimpacts = 40 * 60 * 20;
		VAR_n_save_wave = 1000;
	end

	% BATCH
	BATCH0_wave_damping_scale = VAR_wave_damping_scale;
	BATCH1_h1 = VAR_h1;

	% Saving
	VAR_outputFolder = "RES";

	%% ================================================================

	count0 = length(BATCH0_wave_damping_scale);
	count1 = length(BATCH1_h1);
	threadCount = count0 * count1;

	% Only do one run if using on local
	if isfile(BASE_DIRECTORY + "/ISLOCAL")
		threadCount = 1;
	end

	parfor i = 1:threadCount
		%% Unpack Dispatch Parameters
		idx0 = mod((i - 1), count0) + 1;
		idx1 = floor((i - 1) / count0) + 1;

		%% Set Fluid Parameters

		p = fluid_params();

		%% Set Topography Parameters

		% Domain Size
		domainWidth = VAR_domainWidth;
		p.Lx = domainWidth; % lambdaF
		p.Ly = domainWidth; % lambdaF

		% Grid Spacing
		gridPerWave = VAR_gridPerWave;
		p.hx_desired = 1/gridPerWave;
		p.hy_desired = 1/gridPerWave;
		p.hx         = p.Lx/ceil(p.Lx/p.hx_desired); % lambdaF
		p.hy         = p.Ly/ceil(p.Ly/p.hy_desired); % lambdaF

		p.Nx = p.Lx/p.hx; % dimensionless
		p.Ny = p.Ly/p.hy; % dimensionless

		% Note:
		% round desired grid spacing to whole number of points per lambdaF

		% Time Step (for wave evolution)
		timeStepDivisor = VAR_timeStepDivisor;
		p.dt_desired = min(p.hx,p.hy)^2/timeStepDivisor;    % TF
		p.dt         = 1/ceil(1/p.dt_desired); % TF

		p.nsteps_impact = 1/p.dt; % dimensionless

		% Note:
		% currently set based on spatial resolution, could set independently
		% does not affect drop since instantaneous impacts, only wave evolution
		% round desired time step to have whole number of steps per TF
		% scale with spatial res squared to keep well behaved for high res

		% Topography
		p.topography_type = VAR_topography_type; % options: 'flat', 'square_well', 'circular_well'
		p.damping_type = VAR_damping_type; % options: 'none', 'scaled'
		p.corral_type = VAR_corral_type; % options: 'none', 'rigid', 'spring'

		radius = VAR_R;
		switch p.topography_type
			case "flat"
				p.h1 = VAR_h0_base;
				p.h0 = VAR_h0_base; % m (constant depth)
				p.Rc = radius;      % lambdaF (droplet corral radius)
				p.Dc = p.Rc*2;      % lambdaF (droplet corral diameter)
			case "square_well"
				p.h0 = 1.5*10^(-3); % m (interior depth)
				p.h1 = 1.5*10^(-4); % m (exterior depth)
				p.Lt = 4;           % lambdaF (well width)
			case "circular_well"
				p.h1 = BATCH1_h1(idx1);     % m (interior depth)
				p.h0 = VAR_h0_base + p.h1; % m (exterior depth)
				p.Rc = radius;  % lambdaF (well radius)
				p.Dc = p.Rc*2;  % lambdaF (well diameter)
		end

		switch p.damping_type
			case "scaled"
				p.effective_corral_radius = radius * VAR_effective_corral_radius_scale;
				p.damping_scale = VAR_damping_scale;
		end

		p.wave_damping_scale = BATCH0_wave_damping_scale(idx0);
		p.wave_damping_blend_rate = VAR_wave_damping_blend_rate;

		switch p.corral_type
			case "spring"
				p.effective_corral_radius = radius * VAR_effective_corral_radius_scale;
				p.spring_force_coefficient = VAR_spring_force_coefficient;
		end

		p = top_params(p);

		%% Calculate Faraday Threshold
		if VAR_shouldOverrideThreshold
			fprintf("%s: Overriding threshold.\n", datetime);
			p.GamF = VAR_thresholdGuess;
		else
			thresholdFile = sprintf("%s/threshold_cache/%f_%f_%f_%d_%d_%d_%d_%f.mat", BASE_DIRECTORY, p.h0, p.h1, p.Rc, timeStepDivisor, gridPerWave, domainWidth, p.wave_damping_scale, p.wave_damping_blend_rate);
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

		% Number of Drops
		p.n_drops = VAR_n_drops;

		% Note:
		% only the mass matters since treated as a point for impacts

		% Drop Size
		p.drop_radius  = VAR_r;                                   % m
		p.drop_density = 949;                                     % kg/m^3
		p.drop_mass    = (4/3)*pi*p.drop_radius^3*p.drop_density; % kg

		% Note:
		% effectively controls speed of drop given other parameters
		% optional: can choose a phase to match speed shown in experiments

		p.droplet_collision_type = VAR_droplet_collision_type;
		switch (p.droplet_collision_type)
			case "spring"
				p.droplet_collision_k = VAR_droplet_collision_k;
			case "exp"
				p.droplet_collision_k = VAR_droplet_collision_k;
				p.droplet_collision_n = VAR_droplet_collision_n;
		end

		% Impact Phase
		p.theta     = VAR_theta * pi;

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

		p.c4 = VAR_c4;
		p = drop_params(p);

		% Save Parameters
		p.n_save_wave = VAR_n_save_wave;

		%% Run Simulation

		fprintf("%s: Beginning simulation with time step divisor %d.\n", datetime, timeStepDivisor);
		p = simulate(p);

		%% Output Results

		outputSubfolder = sprintf("RES_N=%d, mem=%.2f, %s R=%.2f h0=%.2f h1=%.2f, theta=%.2f, wds=%d, wdbr=%.2f", p.n_drops, p.mem * 100, p.topography_type, p.Rc, p.h0 * 1000, p.h1 * 1000, p.theta / pi, p.wave_damping_scale, p.wave_damping_blend_rate);
		outputFolder = fullfile(VAR_outputFolder, outputSubfolder);
		if ~isfolder(outputFolder)
			mkdir(outputFolder);
		end
		outputFileName = sprintf("RES_%d.mat", IN_batchIndex);
		saveFilePath = fullfile(outputFolder, outputFileName);
		fprintf("%s: Saving simulation results for %s.\n", datetime, saveFilePath);
		parsave(saveFilePath, p);
	end
end

%% Hack to allow saving inside parfor
function parsave(fname, p)
save(fname, 'p', "-v7.3");
end

function parsave_gamf(fname, GamF)
save(fname, 'GamF', "-v7.3");
end