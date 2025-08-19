function [ui, vi, phi_hat] = drop_impact(xi,yi,ui,vi,phi_hat,eta_hat,p)

for k = 1:p.n_drops

    % Wave Field Gradient at Drop Location
    [Fx,Fy] = compute_slope_eta(eta_hat,xi(k),yi(k),p);

    % Update Drop Speed Due to Instantaneous Impact & For Bounds
    xPos = xi(k);
    yPos = yi(k);
    radius = sqrt(xPos^2 + yPos^2);
    radius_uv = [xPos, yPos] / radius;
    gradient_v = [Fx, Fy];
    velocity_v = [ui(k), vi(k)];

    switch p.damping_type
        case 'scaled'
            % TODO: define p.damping_scale
            % TODO: define p.effective_corral_radius
            if radius > p.effective_corral_radius
                impactScale = p.damping_scale * p.cf_impact;
            else
                impactScale = p.cf_impact;
            end
        otherwise
            impactScale = p.cf_impact;
    end

    % Update Drop Speed Due to Instantaneous Impact
    factor = (p.G/impactScale)*(1-exp(-impactScale)) + exp(-impactScale);
    velocity_v = -gradient_v * factor * velocity_v;

    switch p.corral_type
        case 'rigid'
            if (radius > p.Rc) && (dot(radius_uv, velocity_v) > 0)
                correctedVelocity_v = velocity_v - 2*dot(radius_uv, velocity_v) * radius_uv;

                velocity_v = correctedVelocity_v(1);
                velocity_v = correctedVelocity_v(2);
            end
        case 'spring'
            % TODO: define p.effective_corral_radius
            % TODO: define p.spring_force_coefficient
            if radius > p.effective_corral_radius
                springForce = p.spring_force_coefficient * (radius - p.effective_corral_radius); 
                correctedVelocity_v = velocity_v + spring_force * (-radius_uv);

                velocity_v = correctedVelocity_v(1);
                velocity_v = correctedVelocity_v(2);
            end
    end

    ui(k) = velocity_v(1);
    vi(k) = velocity_v(2);
    
    % Update Velocity Potential Due to Instantaneous Impact
    phi_hat = phi_hat - (p.M*p.G/(p.hx*p.hy))*exp(-p.Kx.*(p.Lx./2+xi(k))-p.Ky.*(p.Ly./2+yi(k)));

end