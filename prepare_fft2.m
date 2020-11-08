function [k_x,k_y,kx,ky,kxx,kyy,k2,k4] = prepare_fft2(Lx,Ly,Nx,Ny)

% k_x = 1i*[0:Nx/2-1 0 -Nx/2+1:-1]*(2*pi/Lx);
% k_y = 1i*[0:Ny/2-1 0 -Ny/2+1:-1]*(2*pi/Ly);

k_x = 1i*[0:Nx/2 -Nx/2+1:-1]*(2*pi/Lx);
k_y = 1i*[0:Ny/2 -Ny/2+1:-1]*(2*pi/Ly);

[kx,  ky ] = meshgrid(k_x,k_y);

k2x = k_x.^2;
k2y = k_y.^2;
[kxx, kyy] = meshgrid(k2x,k2y);

k2 = kxx + kyy;
k4 = k2.^2;

end %endfunction
