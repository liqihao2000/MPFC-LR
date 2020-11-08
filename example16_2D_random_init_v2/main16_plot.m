clear;
clc;
close all;

dirname = 'ex16_2_1_MPFCdata';
% dirname = 'ex16_2_2_MPFCdata';
% dirname = 'ex16_2_3_MPFCdata';

datadir = [dirname,'/data'];
figdir  = [dirname,'_'];

tsave = 10;
nsave = 1000;

Nx = 128; Ny =128;
domain.left   = 0;
domain.right  = 128;
domain.bottom = 0;
domain.top    = 128;
Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;
hx = Lx/Nx;
hy = Ly/Ny;

x  = domain.left   + hx*(0:Nx-1);
y  = domain.bottom + hy*(0:Ny-1);

[xx,yy] = meshgrid(x,y);

kk=1;
% for t= 0:tsave:nsave
for t= [20 100 500 2000]
% for t= [2000]
    figure(kk)
    t
    filename = [datadir '/phi_t=' num2str(t)];
    ss = [filename '.txt'];
    phi = load(ss);
    n   = size(phi(:),1);
    n   = round(n^(1/2));
    phi = reshape(phi,n,n);
    
    mesh(xx,yy,phi,'FaceColor','interp', 'EdgeColor','interp');
    zlim([-1, 1]);
    view(0,90);
    % view(10,13);
    colormap jet;
    axis square;
    axis tight;
    axis off;
    colorbar('Position',[0.845 0.18 0.03 0.66],'Fontsize',15);
%     colorbar off;    

    if t==2000 && (1==strcmp(dirname,'ex16_2_1_MPFCdata'))
        text(62,70,10,'Hexagon','FontSize',30,'FontWeight','bold');
        annotation('arrow',[0.62,0.55],[0.52,0.35],'LineWidth',3)
        text(2,70,10,'Stripes','FontSize',30,'FontWeight','bold');
        annotation('arrow',[0.35,0.65],[0.58,0.80],'LineWidth',3)
    end
    
    if t==2000 && (1==strcmp(dirname,'ex16_2_3_MPFCdata'))
        text(52,30,10,'Hexagon','FontSize',30,'FontWeight','bold');
%       annotation('arrow',[0.62,0.55],[0.52,0.35],'LineWidth',3)
        line([63 56.6 55.6 61.3 67.7 68.8 63]-26,[73.4 70.3 62.9 59 62 69 73.4]-10.2,[1 1 1 1 1 1 1],'Color','k','LineStyle','-','LineWidth',3);
        axes('Position',[0.523,0.583,0.34,0.34]);
        mesh(xx,yy,phi,'FaceColor','interp', 'EdgeColor','interp');
%         mesh(xx(10:50,10:50),yy(10:50,10:50),phi(10:50,10:50),'FaceColor','interp', 'EdgeColor','interp');
        dd = 10.5;
        xlim([37-dd 37+dd]);
        ylim([54-dd 54+dd]);
%         zlim([0,0.4])
        view(0,90);
        
        ax = gca;        
        set(gca,'linewidth',3)
%         ax.XColor = 'c';
%         ax.YColor = 'c';
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        box on        
        axis square;
%         caxis([0 0.4]);
        
        annotation('rectangle',[0.335,0.385,0.107,0.145],'LineWidth',2)
        annotation('arrow',[0.4,0.562],[0.53,0.75],'LineWidth',3)
%         
        line([63 56.6 55.6 61.3 67.7 68.8 63]-26,[73.4 70.3 62.9 59 62 69 73.4]-10.2,[1 1 1 1 1 1 1],'Color','k','LineStyle','-','LineWidth',4);
        
    end
        
    drawnow;
    kk = kk +1;

    
%     figname = [figdir '/phi_t=' num2str(t) '.eps'];
%     print(figname,'-depsc2', '-r120')

    figname = ['../figure_MPFC_SAV/',figdir '_phi_t=' num2str(t) '.jpg'];
    print(figname,'-djpeg', '-r300')
end
