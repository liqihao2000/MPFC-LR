close all;
clear; clc;

% add path
addpath('../','-begin');

N = 64;

% Parameters
para.epsilon = 0.56; 
para.M = 1;
para.beta_bar = 0.1;
para.beta = 1;
para.alpha = 0;
para.sigma = 0.1;
para.name = 'ex19_2_2_MPFCdata';

% Space: Domain and N
L = 50;
domain.xa  = 0;
domain.xb  = L;
domain.ya  = 0;
domain.yb  = L;
domain.za  = 0;
domain.zb  = L; 

% Time: dt T
T = 5000;
t0 = 0;
tsave = 500;

dt_ref = 0.02;

%% option
option.plotflag  = 0;
option.printflag = 1;
option.vtkflag  = 0;
option.saveflag  = 1;
option.savefinal  = 0;
option.energyflag = 1;
option.tol = 1e-14;
option.tolit = 1e-12;
option.maxit = 2000;

pde = ex19_2_MPFCdata(para);

%% Run:

time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
MPFC_3D_CN_SAV(pde,domain,N,N,N,time,option); %% ok

figure(1);
% subplot(1,2,1)
energy=load([pde.name,'e',num2str(pde.epsilon),'M',num2str(pde.M),'b_bar',num2str(pde.beta_bar),'b',num2str(pde.beta),'dt',num2str(time.dt),'_energy.txt']);
% energy = energy(2:40001,:);
% plot(energy(:,1),energy(:,2),'k:',energy(:,1),energy(:,3),'b-','linewidth',1.8);
% fig=plot(energy(:,1),energy(:,3),'b-','linewidth',1.8);
% hold on;
plot(energy(:,1),energy(:,4),'b','markersize',10,'linewidth',1.8);
xlabel('Time','Fontsize',18);ylabel('Energy $E_{cn2}$','Fontsize',18,'interpreter','latex');
set(gca,'FontSize',18);
% set(gca,'YTick',0:100:600);
% ylim([350 550]);
% h = legend('$E$','$\mathcal{E}$','$E_{cn2}$');
% h = legend('$\mathcal{E}$','$E_{cn2}$');
% set(h,'interpreter','latex');
box on;
grid on;
set(gca,'linewidth',1.1)

axes('Position',[0.4,0.45,0.47,0.43]); % subplot
fig2=plot(energy(:,1),energy(:,4),'b','markersize',10,'linewidth',1.8);
annotation('arrow',[0.33,0.20],[0.55,0.48])
xlim([0 300]);
grid on;
set(gca,'linewidth',1.1)
ylim([2500 5500]);
% set(fig2,'YTick',3500:500:6500);
% set(fig2,'xlim',[0 20]);
% axis square

figname = ['../figure_MPFC_SAV/',pde.name,'.eps'];
print(figname,'-depsc2', '-r300')
