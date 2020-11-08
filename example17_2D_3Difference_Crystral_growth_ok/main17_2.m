

close all;
clear; clc;

% add path
addpath('../','-begin');

N = 1024;

% Parameters
para.epsilon = 0.25;
para.M     = 1;
para.beta_bar = 0.1;
para.beta  = 1;
para.alpha = 0;
para.sigma = 0.2;
para.name = 'ex17_2_MPFCdata';

% Space: Domain and N
domain.left   = 0;
domain.right  = 800;
domain.bottom = 0;
domain.top    = 800;

% Time: dt T
T = 3000;
t0 = 0;
tsave = 100;

dt_ref = 0.05;

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

pde = ex17_2_MPFCdata(para);

%% Run:

time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
MPFC_2D_CN_SAV(pde,domain,N,N,time,option); %% ok

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
% set(gca,'XTick',[0:500:3000])
ax = gca;
ax.YAxis.Exponent = 3;%  


% axes('Position',[0.5,0.55,0.35,0.32]); % subplot
% fig2=plot(energy(:,1),energy(:,4),'b','markersize',10,'linewidth',1.8);
% annotation('arrow',[0.45,0.28],[0.7,0.59])
% xlim([0 100]);
% ylim([3500 6500]);
% set(fig2,'YTick',3500:500:6500);
% set(fig2,'xlim',[0 20]);
% axis square

figname = ['../figure_MPFC_SAV/',pde.name,'.eps'];
print(figname,'-depsc2', '-r300')
