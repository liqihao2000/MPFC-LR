%% test order of convergence in space

close all;
clear; clc;

% add path
addpath('../','-begin');

N = 128;

% Parameters
para.epsilon = 0.025;
para.M = 1;
para.beta_bar = 1;
para.beta = 2;
para.alpha = 0;
para.sigma = 0.05;

% Space: Domain and N
domain.left   = 0;
domain.right  = 32;
domain.bottom = 0;
domain.top    = 32;

% Time: dt T
T = 100;
t0 = 0;
tsave = 2;

dt_array = [25 20 10 5 2 1 0.1 0.01 ]';
% dt_array = [25 20 10 5 2 1 ]';

maxIt = length(dt_array);


%% option
option.plotflag  = 0;
option.printflag = 1;
option.vtkflag  = 0;
option.saveflag  = 0;
option.savefinal  = 0;
option.energyflag = 1;
option.tol = 1e-14;
option.tolit = 1e-11;
option.maxit = 2000;

pde = ex14_fig3_MPFCdata(para);

%% Run:
for k = 1:maxIt
    time = struct('T',T,'t0',t0,'dt',dt_array(k),'tsave',tsave);
    MPFC_2D_CN_SAV(pde,domain,N,N,time,option); %% ok
end    

figure(1);
hold on;
lineType = {'>-', 's-','*-','o-','+-' ,'.-' ,'--k','-b'};
for k = 1:maxIt
    energy=load([pde.name,'e',num2str(pde.epsilon),'M',num2str(pde.M),'b_bar',num2str(pde.beta_bar),'b',num2str(pde.beta),'dt',num2str(dt_array(k)),'_energy.txt']);
    plot(energy(:,1),energy(:,4),char(lineType(k)),'LineWidth',1.2);
end
h= legend('$\delta t = 25$','$\delta t = 20$','$\delta t = 10$','$\delta t = 5$','$\delta t = 2$','$\delta t = 1$','$\delta t = 0.1$','$\delta t = 0.01$');
xlabel('Time','Fontsize',18);ylabel('Energy $E_{cn2}$','Fontsize',18,'interpreter','latex');
set(gca,'FontSize',18);
set(gca,'linewidth',1.1)
set(h,'interpreter','latex');
box on;
grid on;
figname = ['../figure_MPFC_SAV/',pde.name,'.eps'];
print(figname,'-depsc2', '-r300')
