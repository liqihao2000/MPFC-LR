%% test of SAV vs IEQ

close all;
clear; clc;

% add path
addpath('../','-begin');

N = 128;

M_array = [5 ];

for kk = 1:length(M_array)
    % Parameters
    para.epsilon = 0.025;
    para.M = M_array(kk);
    para.beta_bar = 1.2;
    para.beta = 0.9;
    para.alpha = 0;
    para.sigma = 2;
    
    % Space: Domain and N
    domain.left   = 0;
    domain.right  = 128;
    domain.bottom = 0;
    domain.top    = 128;
    
    % Time: dt T
    T = 10;
    t0 = 0;
    tsave = 1;
    
    dt_ref = 0.001;
    dt_array = 1./2.^[ 1 2 3 4 5 6 7 8]';
    
    maxIt = length(dt_array);
    
    %% option
    option.plotflag  = 0;
    option.printflag = 0;
    option.vtkflag  = 0;
    option.saveflag  = 0;
    option.savefinal  = 1;
    option.energyflag = 0;
    option.tol = 1e-14;
    option.tolit = 1e-12;
    option.maxit = 2000;
    
    pde = ex03_MPFCdata(para);
    
    %% Run:
    if ~isfield(pde,'exact')
        time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
        MPFC_2D_CN_IEQ(pde,domain,N,N,time,option);
    end
    
    fprintf('T = %.0f\n',T);
    fprintf('  dt   A.n_newton  A.n_cg   CputTime \n');    
    for k = 1:maxIt
        dt = dt_array(k);
        time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
        cput = cputime;
        
        [v2, r1,n_newton,n_cg] = MPFC_2D_CN_IEQ(pde,domain,N,N,time,option);
        
        cput = cputime - cput;
        fprintf('       &  %.2f  &  %.2f  &  %.2f \n',n_newton,n_cg,cput);
    end
    
    %% Compute order of convergence
    error=zeros(maxIt,1);
    order=zeros(maxIt,1);
    if ~isfield(pde,'exact')
        filename=[pde.name,'e',num2str(pde.epsilon),'M',num2str(pde.M),'b_bar',num2str(pde.beta_bar),'b',num2str(pde.beta),'Nx=',num2str(N),'Ny=',num2str(N),'dt=',num2str(dt_ref),'.mat'];
        load(filename,'phi');
        phi_exact = phi;
        clear phi;
    else
        Lx = domain.right - domain.left;
        Ly = domain.top   - domain.bottom;
        hx = Lx/N;
        hy = Ly/N;
        x  = domain.left   + hx*(0:N-1);
        y  = domain.bottom + hy*(0:N-1);
        [xx,yy] = meshgrid(x,y);
        phi_exact = pde.exact(xx,yy,T);
    end
    for k = 1:maxIt
        dt = dt_array(k);
        filenamek=[pde.name,'e',num2str(pde.epsilon),'M',num2str(pde.M),'b_bar',num2str(pde.beta_bar),'b',num2str(pde.beta),'Nx=',num2str(N),'Ny=',num2str(N),'dt=',num2str(dt),'.mat'];
        load(filenamek,'phi');
        error(k,1) = sqrt(sum(sum((phi_exact - phi).^2))*2*pi/N*2*pi/N);   % L2
        clear v;
    end
    order(2:maxIt) = log(error(1:maxIt-1)./error(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));
    
    %% Display error and order
    fprintf('   dt    &   Error_L2\t &  Order \n');
    for k = 1:maxIt
        fprintf('%.4e  %.4e  %.2f \n',dt_array(k),error(k),order(k));
    end
    fprintf('\n')
    
    %% Plot
    hh=loglog(dt_array,error,'*-');
    xlabel('Time step $\delta t$','Interpreter','latex');
    ylabel('$L^2$ error','Interpreter','latex');
    grid on;
    hold on;
    
    %% Save error and order
    name=[pde.name,'e',num2str(pde.epsilon),'M',num2str(pde.M),'b_bar',num2str(pde.beta_bar),'b',num2str(pde.beta),'Nx',num2str(N),'Ny',num2str(N)];
    fileID = fopen([name,'.txt'],'w');
    % fprintf(fileID,'%6s\n','%% Results');
    % fprintf(fileID,'%6s\n','% dt	   &   Error_L2	   &  Order');
    % A = [dt_array error];
    % fprintf(fileID,'%.12f   %.4e   \n',A');
    fprintf(fileID,'%.12f     %.4e      %.2f \n',[dt_array,error,order]');
    fclose(fileID);

end

%% results
% T = 10
%   dt   A.n_newton  A.n_cg   CputTime 
%        &  1.00  &  9.10  &  3.17 
%        &  1.00  &  6.05  &  5.53 
%        &  1.00  &  4.50  &  7.23 
%        &  1.00  &  3.64  &  14.41 
%        &  1.00  &  2.88  &  23.81 
%        &  1.00  &  2.78  &  46.88 
%        &  1.00  &  2.00  &  83.19 
%        &  1.00  &  2.00  &  163.48 
% dt     &   Error_L2	   &  Order 
% 5.0000e-01  1.1749e-01  0.00 
% 2.5000e-01  2.4289e-02  2.27 
% 1.2500e-01  5.9808e-03  2.02 
% 6.2500e-02  1.4902e-03  2.00 
% 3.1250e-02  3.7181e-04  2.00 
% 1.5625e-02  9.2843e-05  2.00 
% 7.8125e-03  2.3196e-05  2.00 
% 3.9063e-03  5.7971e-06  2.00 

