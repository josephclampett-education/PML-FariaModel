function is_stable = faraday_evolve(p,Gam)

% Initial Wave Conditions
pert = 10^(-5);
if isfield(p, "threshold_perturbation_type")
    switch p.threshold_perturbation_type
        case 'bessel'
            eta  = pert.*besselj(0,2*pi.*sqrt(p.xx.^2+p.yy.^2));
        case 'gaussians'
            eta = zeros(size(p.xx));
            num_deltas = 10;
            for i = 1:num_deltas
                sigma = 0.1;

                r = p.Rc * sqrt(rand());
                theta = 2 * pi * rand();
                x = r .* cos(theta);
                y = r .* sin(theta);

                r = sqrt((p.xx - x).^2 + (p.yy - y).^2);
                eta = eta + pert .* exp(-r.^2 / (2 * sigma^2));
            end
    end
else
    eta  = pert.*besselj(0,2*pi.*sqrt(p.xx.^2+p.yy.^2));
end
phi  = zeros(size(p.xx));

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
        fprintf("%s: eta_max %f >> perturbation for Gam %.2f at n = %d. Concluding unstable.\n", datetime, eta_max, Gam, n);
        is_stable = 0;
        break
    elseif eta_max < pert/10
        fprintf("%s: eta_max %f << perturbation for Gam %.2f at n = %d.  Concluding stable.\n", datetime, eta_max, Gam, n);
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