function [xi, yi, ui, vi] = evolve_drops(xi,yi, ui, vi, p)

% Evolve Drop Velocity During Free Flight
ui = ui.*exp(-p.cf_air*p.impact_interval);
vi = vi.*exp(-p.cf_air*p.impact_interval);

% Evolve Drop Position During Free Flight
xi = xi + ui./p.cf_air.*(1-exp(-p.cf_air*p.impact_interval));
yi = yi + vi./p.cf_air.*(1-exp(-p.cf_air*p.impact_interval));
xi = mod(xi+p.Lx/2,p.Lx)-p.Lx/2;
yi = mod(yi+p.Ly/2,p.Ly)-p.Ly/2;

  % Note
  % last two lines force drop to stay in domain, so check to make sure
  % we are getting accurate results if it ever reaches the domain boundary