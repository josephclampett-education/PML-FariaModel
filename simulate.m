function p = simulate(p)

% Set Initial Conditions
t   = p.theta/(4*pi);
phi = p.phi0;
eta = p.eta0; 
xi  = p.xi; 
yi  = p.yi; 
ui  = p.ui; 
vi  = p.vi; 

% Create Variables to Store Data
x_data   = zeros(p.nimpacts,p.n_drops); 
y_data   = zeros(p.nimpacts,p.n_drops); 
t_data   = zeros(p.nimpacts,1);
eta_data = zeros(p.Nx,p.Ny,p.n_save_wave); 

% Fourier Transform of Wave Initial Condition
phi_hat = fft2(phi);               
eta_hat = fft2(eta);

for n = 1:p.nimpacts
    
    % Wave Field at Current Time
    eta = real(ifft2(eta_hat));  

    % Check to Make Sure Not Above Faraday Threshold
    eta_max = max(max(abs(eta)));
    if eta_max > 1
        fprintf("%s: eta_max %f > 1 during simulation.\n", datetime, eta_max);
        % break
    end

    % Store Data
    x_data(n,:)     = xi; 
    y_data(n,:)     = yi; 
    t_data(n)       = t;
    
    % Save wavefield for the last 10 periods
    if n > p.nimpacts - p.n_save_wave
        eta_data(:,:, n - p.nimpacts + p.n_save_wave) = eta;
    end
    
    % Drop Impact
    [ui, vi, phi_hat] = drop_impact(xi,yi, ui, vi, phi_hat, eta_hat, p);
    
    % Evolve Drops Between Impacts
    [xi, yi, ui, vi] = evolve_drops(xi, yi, ui, vi, p);

    % Evolve Wave Between Impacts
    [phi_hat, eta_hat] = evolve_wave(phi_hat, eta_hat, t, p.Gam, p); 
    
    % Update Time
    t = t+p.impact_interval;

end

% Put Data in Output Structure
p.x_data   = x_data; 
p.y_data   = y_data; 
p.t_data   = t_data;
p.eta_data = eta_data; 