function [Fx,Fy] = compute_slope_eta(eta_hat,xi,yi,p)

% Surface Gradient
eta_x_hat = p.Kx.*eta_hat;
eta_y_hat = p.Ky.*eta_hat;

% Shift Fourier Spectrum To Nearest Grid Point
ix = find(p.xx(1,:)>xi,1); 
iy = find(p.yy(:,1)>yi,1);

if isempty(ix); ix=p.Nx; end
if isempty(iy); iy=p.Ny; end

shiftx = p.xx(1,ix)-xi;
shifty = p.yy(iy,1)-yi;

% Compute Gradient at Actual Drop Position
eta_x = real(ifft2(exp(-p.Kx.*shiftx-p.Ky.*shifty).*eta_x_hat));
eta_y = real(ifft2(exp(-p.Kx.*shiftx-p.Ky.*shifty).*eta_y_hat));

Fx = eta_x(iy,ix); 
Fy = eta_y(iy,ix);

  % Note:
  % The shift we do above is undone by multiplying by the exponential
  % before taking the inverse transform. This makes the code more efficient