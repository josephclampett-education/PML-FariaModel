function p = fluid_params()

% Kinematic Viscosity
p.nu = 2*10^(-5); % m^2/s; 

  % Note:
  % only effective value used below

% Effective Kinematic Viscosity
p.nu = 0.865*p.nu; % m^2/s; 

  % Note:
  % set by matching experimental and simulation Faraday threshold
  % choose 0.865nu for 70 Hz, 1.5 mm depth
  % choose 0.803nu for 80 Hz, 6.1 mm depth

% Surface Tension
p.sig = 0.0206; % kg/s^2

% Density
p.rho = 949; % kg/m^3

% Angular Frequency (of bath vibration)
p.omega0 = (2*pi)*70; % 1/s

% Gravity
p.g0 = 9.81; % m/s^2

% Air Viscosity
p.mu_air = 1.8*10^(-5); % kg/(m*s)

% Ideal/Inviscid Dispersion Relation
p.dispEuler = @(k,h)sqrt(tanh(k.*h).*(p.g0.*k+(p.sig/p.rho).*k.^3));