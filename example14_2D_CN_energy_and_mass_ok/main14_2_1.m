%% test order of convergence in space

close all;
clear; clc;

% add path
addpath('../','-begin');

N = 128;

% Parameters
para.epsilon = 0.025;
para.beta_bar = 1;
para.beta = 10;
para.alpha = 0;
para.sigma = 0.1;

% Space: Domain and N
domain.left   = 0;
domain.right  = 32;
domain.bottom = 0;
domain.top    = 32;

% Time: dt T
T = 100;
t0 = 0;
tsave = 2;

dt_ref = 0.01;
dt_array = [ 0.01 0.1 1 2 5 10 ]';

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

pde = ex14_2_1_MPFCdata(para);

%% Run:

time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
MPFC_2D_CN_SAV(pde,domain,N,N,time,option); %% ok

figure(1);
% subplot(2,1,1)
energy=load([pde.name,'e',num2str(pde.epsilon),'M',num2str(pde.M),'b_bar',num2str(pde.beta_bar),'b',num2str(pde.beta),'dt',num2str(dt_ref),'_energy.txt']);
plot(energy(:,1),energy(:,2),'k:',energy(:,1),energy(:,3),'b-','linewidth',1.8);
hold on;
grid on;
plot(energy(1:250:end,1),energy(1:250:end,4),'.r','markersize',18);
ylim([2.6,2.8])
xlabel('Time','Fontsize',18);ylabel('Energy','Fontsize',18);
set(gca,'FontSize',18);
set(gca,'linewidth',1.1)
h = legend('$E$','$\mathcal{E}$','$E_{cn2}$');
set(h,'interpreter','latex');
figname = ['../figure_MPFC_SAV/',pde.name,'.eps'];
print(figname,'-depsc2', '-r300')
% set(gca,'position',[0.15 0.15 .8 .8]);
% subplot(2,1,2)
% mass=load([pde.name,'_mass.txt']);
% plot(mass(:,1),mass(:,2))