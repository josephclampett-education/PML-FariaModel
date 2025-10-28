function [ui, vi, phi_hat] = drop_impact(xi,yi,ui,vi,phi_hat,eta_hat,p)

for k = 1:p.n_drops

    % Wave Field Gradient at Drop Location
    [Fx,Fy] = compute_slope_eta(eta_hat,xi(k),yi(k),p);

    % Update Drop Speed Due to Instantaneous Impact & For Bounds
    position_v = [xi(k), yi(k)];
    radius = norm(position_v);
    position_uv = position_v / radius;
    gradient_v = [Fx, Fy];
    velocity_v = [ui(k), vi(k)];

    switch p.damping_type
        case 'scaled'
            if radius > p.effective_corral_radius
                cf_impact = p.damping_scale * p.cf_impact;
            else
                cf_impact = p.cf_impact;
            end
        otherwise
            cf_impact = p.cf_impact;
    end

    % Update Drop Speed Due to Instantaneous Impact
    velocity_v = (-gradient_v * (p.G/cf_impact)*(1-exp(-cf_impact))) + (exp(-cf_impact) * velocity_v);

    switch p.corral_type
        case 'rigid'
            if (radius > p.Rc) && (dot(position_uv, velocity_v) > 0)
                correctedVelocity_v = velocity_v - 2*dot(position_uv, velocity_v) * position_uv;

                velocity_v = correctedVelocity_v;
            end
        case 'spring'
            % TODO: define p.effective_corral_radius
            % TODO: define p.spring_force_coefficient
            if radius > p.effective_corral_radius
                springForce = p.spring_force_coefficient * (radius - p.effective_corral_radius); 
                correctedVelocity_v = velocity_v + springForce * (-position_uv);

                velocity_v = correctedVelocity_v;
            end
    end

    ui(k) = velocity_v(1);
    vi(k) = velocity_v(2);

    if isfield(p, 'droplet_collision_type') && strcmp(p.droplet_collision_type, 'spring')
    
        distanceThreshold = 2 * p.drop_radius / p.lambdaF;

        for i = 1:p.n_drops-1
            for j = i+1:p.n_drops
                position_i_v = [xi(i), yi(i)];
                velocity_i_v = [ui(i), vi(i)];
                position_j_v = [xi(j), yi(j)];
                velocity_j_v = [ui(j), vi(j)];

                distance_v = position_j_v - position_i_v;
                distance = norm(distance_v);
                distance_uv = distance_v / distance;
    
                if distance < distanceThreshold && distance > 0
                    overlapDistance = distanceThreshold - distance;
    
                    force = p.droplet_collision_k * overlapDistance;
    
                    deltaVelocity_v = force * distance_uv;
    
                    velocity_i_v = velocity_i_v - deltaVelocity_v;
                    velocity_j_v = velocity_j_v + deltaVelocity_v;
                end

                ui(i) = velocity_i_v(1); vi(i) = velocity_i_v(2);
                ui(j) = velocity_j_v(1); vi(j) = velocity_j_v(2);
            end
        end
    end
    
    % Update Velocity Potential Due to Instantaneous Impact
    phi_hat = phi_hat - (p.M*p.G/(p.hx*p.hy))*exp(-p.Kx.*(p.Lx./2+xi(k))-p.Ky.*(p.Ly./2+yi(k)));

end