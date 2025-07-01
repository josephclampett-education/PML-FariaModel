function [ui, vi, phi_hat] = drop_impact(xi,yi,ui,vi,phi_hat,eta_hat,p)

for k = 1:p.n_drops

    % Wave Field Gradient at Drop Location
    [Fx,Fy] = compute_slope_eta(eta_hat,xi(k),yi(k),p);
    
    % Update Drop Speed Due to Instantaneous Impact
    ui(k) = -Fx*(p.G/p.cf_impact)*(1-exp(-p.cf_impact)) + exp(-p.cf_impact)*ui(k);
    vi(k) = -Fy*(p.G/p.cf_impact)*(1-exp(-p.cf_impact)) + exp(-p.cf_impact)*vi(k);
    
    % Update Velocity Potential Due to Instantaneous Impact
    phi_hat = phi_hat - (p.M*p.G/(p.hx*p.hy))*exp(-p.Kx.*(p.Lx./2+xi(k))-p.Ky.*(p.Ly./2+yi(k)));

end