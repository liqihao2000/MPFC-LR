function [phi,r0] = MPFC_3D_First_order_SAV(pde,domain,Nx,Ny,Nz,time,option)
% Solve 2D phase field crystal equaiton
%    (1)    phi_tt + alpha\phi_t = \Delta \mu
%    (2)    \mu = (1+\Delta)^2\phi + (\phi^3 - epsilon\phi)
%    (2)    \mu = (1+\Delta)^2\phi + \beta*phi + \phi*(\phi^2 - epsilon - \beta)
% Qi Li
% 05/23/2019

global dt epsilon M k2 k4 beta_bar beta alpha C0 hx hy hz Lx Ly Lz sigma

if ~exist('option','var'), option = []; end
if ~isfield(option,'tol')
    option.tol = 10^-14;   % default tol
end
if ~isfield(option,'tolit')
    option.tolit = 10^-8;   % default tolit
end
if ~isfield(option,'maxit')
    option.maxit = 2000;   % default maxit
end
if ~isfield(option,'plotflag')
    option.plotflag = 0;   
end
if ~isfield(option,'saveflag')
    option.saveflag = 0;  
end
if ~isfield(option,'savefinal')
    option.savefinal = 0;  
end
if ~isfield(option,'printflag')
    option.printflag = 0;   
end
if ~isfield(option,'vtkflag')
    option.printflag = 0;   
end
if ~isfield(option,'energyflag')
    option.energyflag = 0;   
end
if 1 == option.energyflag
    figname_mass   = [pde.name,'e',num2str(pde.epsilon),'M',num2str(pde.M),'b_bar',num2str(pde.beta_bar),'b',num2str(pde.beta),'dt',num2str(time.dt),'_mass.txt'];
    figname_energy = [pde.name,'e',num2str(pde.epsilon),'M',num2str(pde.M),'b_bar',num2str(pde.beta_bar),'b',num2str(pde.beta),'dt',num2str(time.dt),'_energy.txt'];    
    out1 = fopen(figname_mass,'w');
    out2 = fopen(figname_energy,'w');
end

tol = option.tol;
tolit = option.tolit;
maxit = option.maxit;

%%
T  = time.T;
t  = time.t0;
dt = time.dt;
tsave = time.tsave;

dir_fig  = [pde.name '/fig'];
dir_data = [pde.name '/data'];
epsilon  = pde.epsilon;
M       = pde.M;
alpha    = pde.alpha;
beta_bar = pde.beta_bar;
beta     = pde.beta;
sigma   = pde.sigma; 

Lx = domain.xb - domain.xa;
Ly = domain.yb - domain.ya;
Lz = domain.zb - domain.za;
hx = Lx/Nx;
hy = Ly/Ny;
hz = Ly/Nz;
x  = domain.xa + hx*(0:Nx-1);
y  = domain.ya + hy*(0:Ny-1);
z  = domain.za + hy*(0:Nz-1);

[xx,yy,zz] = meshgrid(x,y,z);
phi0 = pde.init(xx,yy,zz);
psi0 = zeros(size(phi0));
nfigure =1;

%% plot initial value
if 1 == option.saveflag
    if ~exist(dir_data,'dir')
        mkdir(dir_data);
    end
    ss = [dir_data '/phi_t=' num2str(t) '.txt'];
    fid = fopen(ss, 'wt');
    fprintf(fid, '%f\n', phi0(:));
    fclose(fid);
end
if 1 == option.plotflag
    if 1 == option.saveflag
        showsolution_3D_2(nfigure,xx,yy,zz,phi0,t,dir_fig);
        %         write_vtk_grid_values(dir_data,x,y,0,phi0);
    else
        showsolution_3D_2(nfigure,xx,yy,zz,phi0,t);
    end
end

[k_x,k_y,k_z,kx,ky,kz,kxx,kyy,kzz,k2,k4] = prepare_fft3(Lx,Ly,Lz,Nx,Ny,Nz);

nplot = round((T-t)/dt);
nsave = round(tsave/dt);

tstart = tic;

C0 = (epsilon+alpha).^2*Lx*Ly*Lz;
r0 = fun_r_init(phi0);

% Initial energy
if 1 == option.energyflag
    calculate_energy(out1,out2,hx,hy,hz,t,phi0,psi0,r0,C0);
end

for nt = 1:nplot
    t = t+dt;
    
    phi_star = phi0;    
    
    % step 1
    H = fun_H(phi_star);
    
    if isfield(pde,'rhs') && isfield(pde,'exact')
        rhs = dt*pde.rhs(xx,yy,zz,t);
    else
        rhs = 0;
    end
    
    g = get_rhs(phi0,psi0,r0,H)+rhs;
    
    psiA = inv_A(lap_diff(H));
    psiB = inv_A(g);    
    
    gamma = -fftn(H.*psiA);
    gamma = gamma(1,1,1)*hx*hy*hz;
    
    % Step 2      
    Hphi = fftn(H.*psiB);
    Hphi = Hphi(1,1,1)*hx*hy*hz/(1+dt*M*gamma/2);
    
    % Step 3
    phi = dt*M/2*Hphi.*psiA + psiB;     
   

    %% update phi0
    r0 = fun_r(phi,phi0,r0,H); 
    psi0 = fun_psi(phi,phi0);
    phi0 = phi;  
    
    if 1 == option.energyflag
        calculate_energy(out1,out2,hx,hy,hz,t,phi0,psi0,r0,C0);
    end          

    if  0 == mod(nt,nsave)
        if 1 == option.printflag
            timeElapsed = toc(tstart);
            fprintf('epsilon=%.3f,t=%.4f/%.f, dt=%.4f, Nx=%d, Ny=%d, timeElapsed=%f\n',epsilon,t,T,dt,Nx,Ny,timeElapsed);
        end
         
        if 1 == option.saveflag
            ss = [dir_data '/phi_t=' num2str(t) '.txt'];
            fid = fopen(ss, 'wt');
            fprintf(fid, '%f\n', phi(:));
            fclose(fid);
        end
                
        nfigure = nfigure +1;
        if 1 == option.plotflag
            if 1 == option.vtkflag
                write_vtk_grid_values(dir_data,x,y,nt,phi0);
            end
            if 1 == option.saveflag
                showsolution_3D_2(nfigure,xx,yy,zz,phi,t,dir_fig);
            else
                showsolution_3D_2(nfigure,xx,yy,zz,phi,t);
            end
        end
    end
    
end

if 1 == option.savefinal
    name=['phi_e',num2str(epsilon),'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];
    filename=[name '.mat'];
    save(filename,'epsilon','x','y','z','Nx','Ny','Nz','dt','T','phi','domain');
%     showsolution_2D(nfigure,xx,yy,phi,t,'./','.fig');
end
if 1 == option.energyflag
    fclose(out1);
    fclose(out2);
end
end

function r = fun_psi(phi,phi0)
global dt
r = (phi-phi0)/dt;
end

function r = fun_r_init(phi)
global C0 hx hy hz
E1 = fftn(F(phi));
r  = sqrt(E1(1,1,1)*hx*hy*hz + C0);
end

function r = fun_r(phi,phi0,r0,H)
global hx hy hz
Hphi0 = fftn(H.*phi0);
Hphi0 = Hphi0(1,1,1)*hx*hy*hz;
Hphi1 = fftn(H.*phi);
Hphi1 = Hphi1(1,1,1)*hx*hy*hz;
g1 = r0 - 1/2*Hphi0;
r = 1/2*Hphi1+g1;
end

function r = fun_H(phi)
global C0 hx hy hz
E1 = fftn(F(phi));
r = F_derivative(phi)./sqrt(E1(1,1,1)*hx*hy*hz+C0);
end

function r = get_rhs(phi0,psi0,r0,H)
global beta_bar beta M dt hx hy hz sigma Lx Ly Lz
Hphi0 = fftn(H.*phi0);
Hphi0 = Hphi0(1,1,1)*hx*hy*hz;
g1 = r0 - 1/2*Hphi0;
gg = dt/(beta_bar+dt*beta)*fftn((beta_bar+dt*beta)/dt*phi0 + beta_bar*psi0);
r = (beta_bar+dt*beta)/dt*phi0 + beta_bar*psi0 + dt*M*lap_diff(H).*g1 + dt*M*sigma*gg(1,1,1)*hx*hy*hz/(Lx*Ly*Lz);
end

function lap=lap_diff(phi)
global k2
lap=real(ifftn((k2.*fftn(phi))));
end

function r = inv_A(phi)
global dt M k2 alpha beta_bar beta sigma
    r = real(ifftn(fftn(phi)./((beta_bar+dt*beta)/dt-dt*M*k2.*(1+k2).^2-dt*M*alpha.*k2+dt*M*sigma)));
end

function [] = calculate_energy(out1,out2,hx,hy,hz,t,phi,psi,r,C0)
global M epsilon alpha k2 Lx Ly Lz sigma
energy1 = hx*hy*hz*sum(sum(sum( 1/4*phi.^4 - epsilon/2*phi.^2+ 1/2*phi.^2 + lap_diff(phi).*phi  + 1/2*lap_diff(phi).^2  + alpha/2*phi.^2)));

Fpsi = fftn(psi)./k2; Fpsi(1,1,1) = 0; psi_inv_lap = real(ifftn(Fpsi));
ene = -0.5/M*hx*hy*hz*sum(sum(sum(psi_inv_lap.*psi)));

phi_hat = fftn(phi);
psi_OK = phi - phi_hat(1,1,1)*hx*hy*hz/(Lx*Ly*Lz);
Fpsi = fftn(-psi_OK)./k2; Fpsi(1,1,1) = 0; psi_OK = real(ifftn(Fpsi));
psi_hat = fftn(-lap_diff(psi_OK).*psi_OK);
ene2 = 0.5*sigma*psi_hat(1,1,1)*hx*hy*hz;

energy2 = energy1 + ene;
energy3 = hx*hy*hz*sum(sum(sum(1/2*phi.^2 + lap_diff(phi).*phi + 1/2*lap_diff(phi).^2 + alpha/2*phi.^2))) + r.^2 - C0 + ene;
mass   = hx*hy*hz*sum(sum(sum( phi )));
fprintf(out1,'%14.6e  %f \n',t,mass);
fprintf(out2,'%14.6e  %f  %f  %f  \n',t,energy1+ene2,energy2+ene2,energy3+ene2);
end

function r = F_derivative(phi0)
global alpha epsilon
    r = phi0.*(phi0.^2 - epsilon - alpha);
end

function r = F(phi)
global alpha epsilon
    r = 1/4*phi.^4 - epsilon/2*phi.^2 - alpha/2*phi.^2;
end
