%% test order of convergence in space

close all;
clear; clc;

% add path
addpath('../','-begin');

N = 128;

M_array = [5 1000 100000 10000000];

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
    dt_array = 1./2.^[ 1 2 3 4 5 6 7 8 9 10 11]';
%     dt_array = 1./2.^[ 1 2 3 4 5 6 7 ]';
    
    maxIt = length(dt_array);
    
    %% option
    option.plotflag  = 0;
    option.printflag = 1;
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
        MPFC_2D_BDF_SAV(pde,domain,N,N,time,option);
    end
    for k = 1:maxIt
        dt = dt_array(k);
        time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
        v2 = MPFC_2D_BDF_SAV(pde,domain,N,N,time,option);
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
    fprintf('dt     &   Error_L2\t   &  Order \n');
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
% alpha = 0
% dt     &   Error_L2	   &  Order 
% 2.0000e-01  1.6675e+00  0.00 
% 1.0000e-01  1.3904e+00  0.26 
% 5.0000e-02  1.1348e+00  0.29 
% 2.0000e-02  1.3808e-04  9.84 
% 1.0000e-02  3.3635e-05  2.04 
% 5.0000e-03  8.1987e-06  2.04 
% 2.0000e-03  1.2961e-06  2.01 

