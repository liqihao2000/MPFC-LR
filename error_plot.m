% Error plot
clear; clc; close all;
format long;
   
dirname={  'example06_2D_CN_in_Time_exact_ok'
           'example06_2D_CN_in_Time_exact_o_2'
           'example09_2D_BDF_in_Time_exact_ok'
           'example09_2D_BDF_in_Time_exact_o_2' };

epsilon = 0.025;
beta_bar = 1.2;
beta = 0.9;


N=128;

M_array = [5 1000  10000000];
n = 0;

for kkk = 1:length(M_array)
    n = n+1;
    figure(n);
    for kk = 1:size(dirname,1)
        name=[char(dirname(kk)),'/ex03_MPFCdatae',num2str(epsilon),'M',num2str(M_array(kkk)),'b_bar',num2str(beta_bar),'b',num2str(beta),'Nx',num2str(N),'Ny',num2str(N)];
        A = load([name,'.txt']);
        dt_array = A(:,1);
        error(:,kk) = A(:,2);
    end
    loglog(dt_array,1/200*dt_array.^2,'k:','linewidth',2);
    hold on;
    grid on;    
    loglog(dt_array,error(:,1),'v-', 'markersize',10,'linewidth',2);
    loglog(dt_array,error(:,2),'X-', 'markersize',10,'linewidth',2);
    loglog(dt_array,error(:,3),'s-', 'markersize',10,'linewidth',2);
    loglog(dt_array,error(:,4),'*-', 'markersize',10,'linewidth',2);
    
    ylim([10e-10 10e1])
    legend({'$\mathcal{O}(\delta t^2)$', 'SAV-CN', 'S-SAV-CN','SAV-BDF2', 'S-SAV-BDF2'},'Interpreter','latex','Location','southeast','Fontsize',12);
    xlabel('Time step $\delta t$','Interpreter','latex','Fontsize',14);
    ylabel('$L^2$ error','Interpreter','latex','Fontsize',14);
    set(gca,'FontSize',14);
    set(gca,'linewidth',1.1)
    set (gcf,'position',[300,100,630,478] )  
    set(gca,'YTick',[10e-10,10e-8,10e-6,10e-4,10e-2,10e-0])
    set(gca, 'YMinorGrid','off');
    set(gca, 'YMinorTick','off');
    
%     if n == 1     
%         axes('position', [0.3 0.5 0.5 0.5])
%         grid on;
%         loglog(dt_array,1/200*dt_array.^2,'k:','linewidth',2);
%         hold on;
%         grid on;
%         loglog(dt_array,error(:,1),'v-', 'markersize',10,'linewidth',2);
%         loglog(dt_array,error(:,2),'X-', 'markersize',10,'linewidth',2);
%         loglog(dt_array,error(:,3),'s-', 'markersize',10,'linewidth',2);
%         loglog(dt_array,error(:,4),'*-', 'markersize',10,'linewidth',2);
%         xlim([10e-4, 10e-3 ])
% 
%     end    
    
     
    if n == 2
        annotation('arrow',[0.556,0.556],[0.295,0.375],'linewidth',1.5);
        text(2^(-6)/2,1e-7,'$\delta t = 2^{-6}$','FontSize',13,'Interpreter','latex');
        %             annotation('arrow',[0.555,0.51],[0.56,0.55],'linewidth',1.5);
        %             text(2^(-7)/5,3e-4,'$\delta t = 2^{-7}/10$','FontSize',13,'Interpreter','latex');
        annotation('arrow',[0.27,0.365],[0.43,0.37],'linewidth',1.5);
        text(2^(-7)/30,30e-6,'$\delta t = 2^{-9}$','FontSize',13,'Interpreter','latex');
%         figname = ['../../papers/paper3_2_SAV_MPFC_20190117/0_Submission/figure_MPFC_SAV/error_',num2str(n),'.eps'];
%         print(figname,'-depsc2', '-r600')
    end
    
%     if n == 3
%         annotation('arrow',[0.50,0.497],[0.42,0.49],'linewidth',1.5);
%         text(2^(-11),3e-6/3,'$\delta t = 2^{-7}/10$','FontSize',13,'Interpreter','latex');
%     end
    
%     if n == 5
%         annotation('arrow',[0.32,0.37],[0.81,0.72],'linewidth',1.5);
%         text(2^(-11)/10,4,'$\delta t = 2^{-9}/10$','FontSize',13,'Interpreter','latex');
%         annotation('arrow',[0.495,0.45],[0.53,0.58],'linewidth',1.5);
%         text(2^(-7)/10,1e-4,'$\delta t = 2^{-8}/10$','FontSize',13,'Interpreter','latex');
%         annotation('arrow',[0.73,0.69],[0.56,0.54],'linewidth',1.5);
%         text(2^(-6)/1.2,3e-4,'$\delta t = 2^{-4}/10$','FontSize',13,'Interpreter','latex');
%         annotation('arrow',[0.27,0.31],[0.376,0.305],'linewidth',1.5);
%         text(2^(-13)/10,1e-6,'$\delta t = 2^{-10}/10$','FontSize',13,'Interpreter','latex');
%     end
%     
%     if n == 6
%         annotation('arrow',[0.65,0.57],[0.52,0.50],'linewidth',1.5);
%         text(2^(-6)/3,3e-4/4,'$\delta t = 2^{-6}/10$','FontSize',13,'Interpreter','latex');
%     end

    figname = ['./figure_MPFC_SAV/error_',num2str(n),'.eps'];
    print(figname,'-depsc2', '-r600')
end

