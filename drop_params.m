function p = drop_params(p)

% Dimensionless Drop Mass
p.M = p.drop_mass./(p.rho*p.xF^3);

% Air Drag Coefficient
p.cf_air = 6*pi*p.drop_radius*p.mu_air*p.TF/p.drop_mass;

% Stokes (Fluid) Drag Coefficient
p.c4 = 0.17;
p.cf_impact = p.c4*sqrt(p.rho*p.drop_radius/p.sig)*p.TF*p.g0;

  % Note:
  % c4 is the coefficient of restitution (see Molacek 2013)

% Dimensionless Drop Time Between Impacts
p.impact_interval = 1;

  % Note:
  % once per Faraday period, change if not (2,1) mode

%%% Strobe Model Parameters %%%

% Drag Coefficient
p.strobe.D = p.cf_air + p.cf_impact;

