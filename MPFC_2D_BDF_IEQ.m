function [phi,U0,n_newton,n_cg] = MPFC_2D_BDF_IEQ(pde,domain,Nx,Ny,time,option)
% Solve 2D phase field crystal equaiton
%    (1)    phi_tt + alpha\phi_t = \Delta \mu
%    (2)    \mu = (1+\Delta)^2\phi + (\phi^3 - epsilon\phi)
%    (2)    \mu = (1+\Delta)^2\phi + \beta*phi + \phi*(\phi^2 - epsilon - \beta)
% Qi Li
% 01/16/2019
global dt epsilon M k2 k4 beta_bar beta alpha C0 B hx hy Lx Ly sigma

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
    out1 = fopen(figname_mass,'a');
    out2 = fopen(figname_energy,'a');
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

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;
hx = Lx/Nx;
hy = Ly/Ny;
% x  = domain.left   + hx*(1:Nx);
% y  = domain.bottom + hy*(1:Ny);
x  = domain.left   + hx*(0:Nx-1);
y  = domain.bottom + hy*(0:Ny-1);

[xx,yy] = meshgrid(x,y);
phi0 = pde.init(xx,yy);
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
        showsolution_2D(nfigure,xx,yy,phi0,t,dir_fig);
        %         write_vtk_grid_values(dir_data,x,y,0,phi0);
    else
        showsolution_2D(nfigure,xx,yy,phi0,t);
    end
end

[k_x,k_y,kx,ky,kxx,kyy,k2,k4] = prepare_fft2(Lx,Ly,Nx,Ny);

nplot = round((T-t)/dt);
nsave = round(tsave/dt);

tstart = tic;

% C0 = (epsilon+alpha).^2*Lx*Ly;
B = (epsilon+alpha).^2;
U0 = fun_U_init(phi0);

%% initialization phi 1
time1 = time; time1.T = dt;
[phi1, U1,n_newton,n_cg] = MPFC_2D_First_order_IEQ(pde,domain,Nx,Ny,time1,option);
[phi1, U1] = MPFC_2D_CN_IEQ_PreCor(pde,domain,Nx,Ny,time1,phi1,option);

psi1 = 2*(phi1 - phi0)/dt-psi0; 

t = t+dt;

for nt = 2:nplot
    t = t+dt;
    
    phi_star = 2*phi1 - phi0;
    
    % step 1
    H = fun_H(phi_star);
    
    if isfield(pde,'rhs') && isfield(pde,'exact')
        rhs = dt*pde.rhs(xx,yy,t);
    else
        rhs = 0;
    end
    
    g = get_rhs(phi1,phi0,psi1,psi0,U1,U0,H)+rhs;  
    
    x0        = phi0(:);
    lhs       = @get_lhs;
    inv       = @inv_M;
    [phi,flag,relres,iter]= pcg(@(x)lhs(x,H),g(:),tolit,maxit,@(y)inv(y,H),[],x0);
    phi = reshape(phi,Nx,Ny); 
    
    n_cg = n_cg + iter;    
    
    
%% update phi0  
    Uold = U1;
    U1 = fun_U(phi,phi1,phi0,U1,U0,H);
    U0 = Uold;  

    psi0 = psi1;
    psi1 = fun_psi(phi,phi1,phi0);
    
    phi0 = phi1;      
    phi1 = phi;     
    
    if 1 == option.energyflag
        calculate_energy(out1,out2,hx,hy,t,phi0,phi1,psi0,psi1,U0,U1,C0)
    end

    if  0 == mod(nt,nsave)
        if 1 == option.printflag
            timeElapsed = toc(tstart);
            fprintf('epsilon=%.3f,t=%.4f/%.f, dt=%.4e, Nx=%d, Ny=%d, timeElapsed=%f\n',epsilon,t,T,dt,Nx,Ny,timeElapsed);
        end
         
        nfigure = nfigure +1;
        if 1 == option.plotflag
            if 1 == option.vtkflag
                write_vtk_grid_values(dir_data,x,y,nt,phi0);
            end
            if 1 == option.saveflag
                showsolution_2D(nfigure,xx,yy,phi,t,dir_fig);
            else
                showsolution_2D(nfigure,xx,yy,phi,t);
            end
        end
    end
    
end

if 1 == option.savefinal
    name= [pde.name,'e',num2str(pde.epsilon),'M',num2str(pde.M),'b_bar',num2str(pde.beta_bar),'b',num2str(pde.beta),'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];
    filename=[name '.mat'];
    save(filename,'epsilon','M','x','y','hx','hy','Nx','Ny','dt','T','phi','domain');
%     showsolution_2D(nfigure,xx,yy,phi,t,'./','.fig');
end
if 1 == option.energyflag
    fclose(out1);
    fclose(out2);
end
n_cg     = n_cg/nplot;
end

function r = fun_psi(phi,phi1,phi0)
global dt
r = (3*phi-4*phi1+phi0)/(2*dt);
end

function r = fun_U_init(phi)
global B
% E1 = fft2(F(phi));
r  = sqrt(F(phi) + B);
end

function r = fun_U(phi,phi1,phi0,U1,U0,H)
g1 = (4*U1-U0)/3 - 1/2*H.*(4*phi1-phi0)/3;
r = 1/2*H.*phi + g1;
end

function r = fun_H(phi)
global B
% E1 = fft2(F(phi));
r = F_derivative(phi)./sqrt(F(phi)+B);
end

function r = get_rhs(phi1,phi0,psi1,psi0,U1,U0,H)
global beta_bar beta dt M hx hy sigma Lx Ly
coe = (3*beta_bar/2+dt*beta)*3/(2*dt)+dt*M*sigma;
g1 = (4*U1-U0)/3 - 1/2*H.*(4*phi1-phi0)/3;
g2 = -(4*phi1-phi0)/(2*dt);
g_hat = -(3/2*beta_bar+beta*dt)*g2 + beta_bar*(4*psi1-psi0)/2;
gg = 2*dt/((3/2*beta_bar+dt*beta)*3)*fft2(g_hat);
g3 = g_hat+dt*M*sigma*gg(1,1)*hx*hy/(Lx*Ly);
g4 = H.*g1;
r  = g3 + dt*M*lap_diff(g4);
end

function r = get_lhs(phi,H)
global dt M beta_bar beta sigma
n   = size(phi(:),1);
n   = n^(1/2);
phi = reshape(phi,n,n);
coe = (3*beta_bar/2+dt*beta)*3/(2*dt)+dt*M*sigma;
P     = (phi+2*lap_diff(phi)+lap_diff(lap_diff(phi))) + 1/2*H.^2.*phi;
r = coe*phi-dt*M*lap_diff(P);
r = r(:);
end

function r = inv_M(phi,H)
global k2 dt beta_bar beta M sigma
n   = size(phi(:),1);
n   = n^(1/2);
phi = reshape(phi,n,n);
coe = (3*beta_bar/2+dt*beta)*3/(2*dt)+dt*M*sigma;
aver  = mean(1/2*H(:).*H(:));
Matirx = coe - dt*M*k2.*((1+k2).^2+aver);
r = real(ifft2(Matirx.\fft2(phi)));
r = r(:);
end

function lap=lap_diff(phi)
global k2
lap=real(ifft2((k2.*fft2(phi))));
end

function [] = calculate_energy(out1,out2,hx,hy,t,phi1,phi,psi1,psi,r1,r,C0)
global M epsilon beta_bar k2 Lx Ly sigma
energy1 = hx*hy*sum(sum( 1/4*phi.^4 - epsilon/2*phi.^2+ 1/2*phi.^2 + lap_diff(phi).*phi  + 1/2*lap_diff(phi).^2 ));

Fpsi = fft2(psi)./k2; Fpsi(1,1) = 0; psi_inv_lap = real(ifft2(Fpsi));
eneA  = -0.5/M*beta_bar*hx*hy*sum(sum(psi_inv_lap.*psi));

psi_B = 2*psi-psi1;
FpsiB = fft2(psi_B)./k2; FpsiB(1,1) = 0; psi_inv_lapB = real(ifft2(FpsiB));
eneB = -0.5/M*beta_bar*hx*hy*sum(sum(psi_inv_lapB.*psi_B ));
ene = 0.5*(eneA+eneB);

phi_hat = fft2(phi);
psi_OK = phi - phi_hat(1,1)*hx*hy/(Lx*Ly);
Fpsi = fft2(-psi_OK)./k2; Fpsi(1,1) = 0; psi_OK = real(ifft2(Fpsi));
psi_hat = fft2(-lap_diff(psi_OK).*psi_OK);
ene2A = 0.5*sigma*psi_hat(1,1)*hx*hy;

phi1_hat = fft2(phi1);
psi_OKB = 2*(phi - phi_hat(1,1)*hx*hy/(Lx*Ly)) -(phi1 - phi1_hat(1,1)*hx*hy/(Lx*Ly)) ;
Fpsi = fft2(-psi_OKB)./k2; Fpsi(1,1) = 0; psi_OKB = real(ifft2(Fpsi));
psi_hat = fft2(-lap_diff(psi_OKB).*psi_OKB);
ene2B = 0.5*sigma*psi_hat(1,1)*hx*hy;
ene2 = 0.5*(ene2A+ene2B);

energy2 = energy1 + ene;
energy3 = hx*hy*sum(sum(    + 0.25*( phi.^2 + 2*lap_diff(phi).*phi  + lap_diff(phi).^2 ...
                         +(2*phi-phi1).^2 + 2*lap_diff(2*phi-phi1).*(2*phi-phi1)  + lap_diff(2*phi-phi1).^2)))...
                         + 0.5*(r.^2 + (2*r-r1).^2) - C0 + ene;                     
mass    = hx*hy*sum(sum( phi ));
fprintf(out1,'%14.6e  %.8f \n',t,mass);
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
