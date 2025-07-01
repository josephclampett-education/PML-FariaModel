function is_stable = faraday_evolve(p,Gam)

% Initial Wave Conditions
pert = 10^(-5);
eta  = pert.*besselj(0,2*pi.*sqrt(p.xx.^2+p.yy.^2));
phi  = zeros(size(p.xx));

  % Note:
  % Could use gaussian instead of bessel perturbation

% Fourier Transform of Initial Conditions
phi_hat = fft2(phi);               
eta_hat = fft2(eta);

% Setup Vector to Track Square Root of Energy
t_vec       = 1:p.nF;
sqrt_energy = zeros(1,p.nF); % ~ exp(a.*tvec), a = 0 at threshold

% Evolve Wave for p.nF Faraday Periods
for n = 1:p.nF

    % Evolve wave for a single Faraday period
    [phi_hat, eta_hat] = evolve_wave(phi_hat,eta_hat,t_vec(n)-1,Gam,p);

    % Check for obvious stability/instability
    eta     = real(ifft2(eta_hat));
    eta_max = max(max(abs(eta)));
    
    if eta_max > 10*pert
        disp('Probably above threshold. Stopping simulation...');
        is_stable = 0;
        break
    elseif eta_max < pert/10
        disp('Probably below threshold. Stopping simulation...');
        is_stable = 1;
        break
    end
    
    % Update square root of energy
    sqrt_energy(n) = sqrt(sum(abs(eta_hat(:)).^2));

    % Check stability if made it to end of for loop
    if n == p.nF

        % The log of the energy should be linear for late times
        log_energy = log(sqrt_energy);

        % Fit to line
        c = polyfit(t_vec(end-29:end),log_energy(end-29:end),1);

        if c(1) < 0
            is_stable = 1;
        else
            is_stable = 0;
        end

    end

end