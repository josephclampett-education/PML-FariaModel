function p = top_params(p)

% Define Spatial Grid
p.x = p.hx*(0:p.Nx-1)-p.Lx/2;
p.y = p.hy*(0:p.Ny-1)-p.Ly/2;

[p.xx,p.yy] = meshgrid(p.x,p.y);

% Define Topography
p.h_top_grid = zeros(p.Ny,p.Nx);

switch p.type
    case 'flat'
        p.h_top_grid = p.h0;
    case 'square_well'
        for i = 1:length(p.x)
            for j = 1:length(p.y)
                xi = p.x(i);
                yj = p.y(j);
                if (abs(xi) <= p.Lt/2) && (abs(yj) <= p.Lt/2)
                    p.h_top_grid(j,i) = p.h0;
                else
                    p.h_top_grid(j,i) = p.h1;
                end
            end
        end
    case 'circular_well'
        for i = 1:length(p.x)
            for j = 1:length(p.y)
                xi = p.x(i);
                yj = p.y(j);
                if sqrt(xi^2+yj^2) <= p.Rc
                    p.h_top_grid(j,i) = p.h0;
                else
                    p.h_top_grid(j,i) = p.h1;
                end
            end
        end
end

% Mean Depth
p.h_mean = mean(mean(p.h_top_grid));

% (Approximate) Wave Number at Each Point
p.wi = p.dispEuler;               % ideal/inviscid dispersion
p.wg = @(k,h)p.g0.*k.*tanh(k.*h); % gravity wave dispersion
p.gm = @(k)2*p.nu.*k.^2;          % dissipation rate

p.Gamma_neutral = @(k,h)sqrt((4./p.wg(k,h).^2).* ... % need to find min (Eq 2.8 Faria)
    ((p.wi(k,h).^2+p.gm(k).^2-(p.omega0/2)^2).^2+(p.omega0.*p.gm(k)).^2));

p.kf_deep        = zeros(size(p.h_top_grid));
p.Gamma_max_deep = zeros(size(p.h_top_grid));

options = optimset('Display','off'); % suppress print statements for min search

for i = 1:length(p.h_top_grid(:,1))
    for j = 1:length(p.h_top_grid(1,:))
        hij = p.h_top_grid(i,j);
        [p.kf_deep(i,j),p.Gamma_max_deep(i,j)] = fminsearch(@(k)p.Gamma_neutral(k,hij),1300,options);
    end
end

[p.kf_mean,p.Gamma_max_mean] = fminsearch(@(k)p.Gamma_neutral(k,p.h_mean),1300,options);

% Effective Faraday Wavelength
p.lambdaF = (2*pi)/p.kf_mean;

% Effective Depth
p.d_deep = tanh(p.kf_deep.*p.h_top_grid)./p.kf_deep;

% Scales for Non-dimensionalization
p.xF = p.lambdaF;           % Faraday wavelength
p.TF = (2*pi)/(p.omega0/2); % Faraday period

% Reynolds Number
p.Reynolds = p.xF.^2/(p.TF*p.nu);
p.nu0      = 1/p.Reynolds;

% Bond Number
p.Bo = p.sig*p.TF^2/(p.rho*p.xF^3);

% Dimensionless Gravity
p.G = p.g0*p.TF^2/p.xF;
p.g = @(t,Gam)p.G.*(1-Gam.*cos(4*pi.*t));

% Dimensionless Depth
p.h0_deep = p.h_top_grid./p.xF;        p.h = p.h0_deep;
p.d0_deep = p.d_deep./p.xF;            p.d = p.d0_deep;        
p.a0_deep = zeros(size(p.h_top_grid)); p.a = p.a0_deep;

% Dimensionless Faraday Wave Number
p.kf0_deep = p.kf_deep.*p.xF;

% Wave Number Grid in Fourier Space
p.kx = 2*pi*1i/p.Lx*[0:p.Nx/2-1 0 -p.Nx/2+1:-1];
p.ky = 2*pi*1i/p.Ly*[0:p.Ny/2-1 0 -p.Ny/2+1:-1];

[p.Kx,p.Ky] = meshgrid(p.kx,p.ky);

% (Why?) Dissipation Operator in Fourier Space
p.K2              = p.Kx.^2 + p.Ky.^2;
p.abs_K           = sqrt(-p.K2);
p.dissMatrix      = exp(2*p.nu0*p.dt*p.K2);      
p.dissMatrix_half = exp(2*p.nu0*(p.dt/2)*p.K2);

% Misc. Definitions Needed to Calculate eta_t
p.shift1 = mod(-(1:p.Nx)+1,p.Nx)+1;
p.shift2 = mod(-(1:p.Ny)+1,p.Ny)+1;
p.KxiKy  = p.Kx+1i*p.Ky;
p.KxmiKy = p.Kx-1i*p.Ky;

% Misc. Definitions Needed to Calculate phi_t
p.kx_deriv              = 2*pi*1i/p.Lx*[0:(p.Nx/2-1) (-p.Nx/2):-1];
p.ky_deriv              = 2*pi*1i/p.Ly*[0:(p.Ny/2-1) (-p.Ny/2):-1];
[p.Kx_deriv,p.Ky_deriv] = meshgrid(p.kx_deriv,p.ky_deriv);
p.K2_deriv              = p.Kx_deriv.^2 + p.Ky_deriv.^2;

% Dissipation Operator in Fourier Space for Wave Evolution
p.D = exp(p.nu0*p.K2_deriv*p.dt);
