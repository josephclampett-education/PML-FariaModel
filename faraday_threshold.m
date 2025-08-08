function p = faraday_threshold(p, guess)

% Starting Guess for Faraday Threshold
Gam = guess; % units of gravitational acceleration g

% Find Upper and Lower Bound of Faraday Threshold
p.nF      = 1000;                  % number of Faraday periods to simulate over
is_stable = faraday_evolve(p,Gam); % returns true if stable, false if unstable

if is_stable
    Gam_lb = Gam;
    Gam_ub = Gam + 0.1;
    while is_stable == true
        is_stable = faraday_evolve(p,Gam_ub);
        if is_stable
            Gam_lb = Gam_lb + 0.1;
            Gam_ub = Gam_ub + 0.1;
        end
    end
else
    Gam_lb = Gam - 0.1;
    Gam_ub = Gam;
    while is_stable == false
        is_stable = faraday_evolve(p,Gam_lb);
        if ~is_stable
            Gam_lb = Gam_lb - 0.1;
            Gam_ub = Gam_ub - 0.1;
        end
    end
end

% Use Bisection to Find Faraday Threshold
for i = 1:10
    Gam       = (Gam_lb+Gam_ub)/2;
    is_stable = faraday_evolve(p,Gam);
    if is_stable
        Gam_lb = Gam;
    else
        Gam_ub = Gam;
    end
end

p.GamF = (Gam_lb+Gam_ub)/2;