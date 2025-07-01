function [phi_hat, eta_hat] = evolve_wave(phi_hat, eta_hat, t_in, Gam, p)

t = t_in;

for n=1:p.nsteps_impact  
    
    % stage 1
    [rhs1_1, rhs2_1] = compute_rhs_full_IF(phi_hat, eta_hat, t, Gam, p);
    rhs1_1 = p.D .* rhs1_1; rhs2_1 = p.D .* rhs2_1;
    
    % stage 2
    [rhs1_2, rhs2_2] = compute_rhs_full_IF(p.D.*phi_hat + p.dt/2.*rhs1_1, ...
        p.D.*eta_hat + p.dt/2.*rhs2_1, t+p.dt/2, Gam, p); 
    
    % stage 3
    [rhs1_3, rhs2_3] = compute_rhs_full_IF(p.D.*phi_hat + p.dt/2.*rhs1_2, ...
        p.D.*eta_hat + p.dt/2.*rhs2_2, t+p.dt/2, Gam, p); 
    
    % stage 4
    [rhs1_4, rhs2_4] = compute_rhs_full_IF(p.D .* (p.D .* phi_hat + p.dt.*rhs1_3), ...
        p.D .* (p.D .* eta_hat + p.dt.*rhs2_3), t+p.dt, Gam, p);

    % RK step
    phi_hat = p.D.^2 .* phi_hat + p.dt/6 * (p.D.*rhs1_1 + 2*p.D.*rhs1_2 + 2*p.D.*rhs1_3 + rhs1_4);
    eta_hat = p.D.^2 .* eta_hat + p.dt/6 * (p.D.*rhs2_1 + 2*p.D.*rhs2_2 + 2*p.D.*rhs2_3 + rhs2_4);

    t = t+p.dt;

end 

end

function [rhs1, rhs2] = compute_rhs_full_IF(phi_hat,eta_hat,t,Gam,p)

    rhs1 = -p.g(t,Gam).*eta_hat + p.Bo*p.K2_deriv.*eta_hat;
    rhs2 = DtN(phi_hat,p);   

end

function [rhs2] =  DtN(phi_hat,p)

% Approximation to \phi_z (i.e. Dirichelt-to-Neumann operator)

w    = p.d.*ifft2(p.KxiKy.*phi_hat);
A    = fft2(w); 
As   = conj(A(p.shift1,p.shift2));
rhs2 = -(p.KxmiKy.*A/2 + p.KxiKy.*As/2);

end