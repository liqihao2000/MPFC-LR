function [k_x,k_y,k_z,kx,ky,kz,kxx,kyy,kzz,k2,k4] = prepare_fft3(Lx,Ly,Lz,Nx,Ny,Nz)

%---

k_x = 1i*[0:Nx/2 -Nx/2+1:-1]*(2*pi/Lx);
k_y = 1i*[0:Ny/2 -Ny/2+1:-1]*(2*pi/Ly);
k_z = 1i*[0:Nz/2 -Nz/2+1:-1]*(2*pi/Lz);

% k_x = 1i*[0:Nx/2-1 0 -Ny/2+1:-1]*(2*pi/Lx);
% k_y = 1i*[0:Ny/2-1 0 -Ny/2+1:-1]*(2*pi/Ly);
% k_z = 1i*[0:Nz/2-1 0 -Nz/2+1:-1]*(2*pi/Lz);

[kx, ky, kz] = meshgrid(k_x,k_y,k_z);

k2x = k_x.^2;
k2y = k_y.^2;
k2z = k_z.^2;

[kxx, kyy, kzz] = meshgrid(k2x,k2y,k2z);

k2 = kxx + kyy + kzz;
k4 = k2.^2;

end %endfunction