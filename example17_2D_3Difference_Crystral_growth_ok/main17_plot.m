clear;
clc;
close all;
dirname = 'ex17_2_MPFCdata';

datadir = [dirname,'/data'];
figdir  = [dirname,'_'];

tsave = 100;
nsave = 1000;

Nx = 1024; Ny =1024;
domain.left   = 0; 
domain.right  = 800;
domain.bottom = 0;
domain.top    = 800;
Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;
hx = Lx/Nx;
hy = Ly/Ny;

x  = domain.left   + hx*(0:Nx-1);
y  = domain.bottom + hy*(0:Ny-1);

[xx,yy] = meshgrid(x,y);

kk=1;
% for t= 0:tsave:nsave
% for t= [ 3000 ]
for t= [ 0 800 1600 3000]
    figure(kk);
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

  if 3000 == t
        axes('Position',[0.523,0.583,0.34,0.34]);
        mesh(xx,yy,phi,'FaceColor','interp', 'EdgeColor','interp');
%         mesh(xx(10:50,10:50),yy(10:50,10:50),phi(10:50,10:50),'FaceColor','interp', 'EdgeColor','interp');
        xlim([37 117]);
        ylim([400 480]);
%         zlim([0,0.4])
        view(0,90);
        
        ax = gca;        
        set(gca,'linewidth',2)
%         ax.XColor = 'c';
%         ax.YColor = 'c';
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        box on        
        axis square;
%         caxis([0 0.4]);
        
        annotation('rectangle',[0.24,0.515,0.06,0.08],'LineWidth',0.8)
        annotation('arrow',[0.3,0.565],[0.554,0.73],'LineWidth',1)
  elseif 1600 == t
        axes('Position',[0.523,0.583,0.34,0.34]);
        mesh(xx,yy,phi,'FaceColor','interp', 'EdgeColor','interp');
%         mesh(xx(10:50,10:50),yy(10:50,10:50),phi(10:50,10:50),'FaceColor','interp', 'EdgeColor','interp');
        xlim([28 108]);
        ylim([410 490]);
%         zlim([0,0.4])
        view(0,90);
        
        ax = gca;        
        set(gca,'linewidth',2)
%         ax.XColor = 'c';
%         ax.YColor = 'c';
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        box on        
        axis square;
%         caxis([0 0.4]);
        
        annotation('rectangle',[0.23,0.53,0.06,0.08],'LineWidth',0.8)
        annotation('arrow',[0.29,0.565],[0.57,0.72],'LineWidth',1)
  elseif 800 == t
        axes('Position',[0.603,0.673,0.25,0.25]);
        mesh(xx,yy,phi,'FaceColor','interp', 'EdgeColor','interp');
%         mesh(xx(10:50,10:50),yy(10:50,10:50),phi(10:50,10:50),'FaceColor','interp', 'EdgeColor','interp');
        xlim([340-30 340+30]);
        ylim([400-30 400+30]);
%         zlim([0,0.4])
        view(0,90);
        box on        
        axis square;
        ax = gca;        
        set(gca,'linewidth',2)
%         ax.XColor = 'c';
%         ax.YColor = 'c';
        set(gca,'XTick',[])
        set(gca,'YTick',[])
%         caxis([0 0.4]);
        annotation('rectangle',[0.45,0.475,0.05,0.07],'LineWidth',0.8)
        annotation('arrow',[0.475,0.632],[0.545,0.74],'LineWidth',1)        

        %%(2)
        axes('Position',[0.182,0.673,0.25,0.25]);
        mesh(xx,yy,phi,'FaceColor','interp', 'EdgeColor','interp');
%         mesh(xx(10:0,10:50),yy(10:50,10:50),phi(10:50,10:50),'FaceColor','interp', 'EdgeColor','interp');
        xlim([190-30 190+30]);
        ylim([205-30 205+30]);
%         zlim([0,0.4])
        view(0,90);
        box on        
        axis square;
        
        ax = gca;        
        set(gca,'linewidth',2)
%         ax.XColor = 'c';
%         ax.YColor = 'c';
        set(gca,'XTick',[])
        set(gca,'YTick',[])
%         caxis([0 0.4]);
        
        annotation('rectangle',[0.31,0.275,0.05,0.07],'LineWidth',0.8)
        annotation('arrow',[0.335,0.30],[0.345,0.67],'LineWidth',1)
  
        %%(3)
        axes('Position',[0.602,0.110,0.25,0.25]);
        mesh(xx,yy,phi,'FaceColor','interp', 'EdgeColor','interp');
%         mesh(xx(10:0,10:50),yy(10:50,10:50),phi(10:50,10:50),'FaceColor','interp', 'EdgeColor','interp');
        xlim([600-30 600+30]);
        ylim([340-30 340+30]);
%         zlim([0,0.4])
        view(0,90);
        box on        
        axis square;
        
        ax = gca;        
        set(gca,'linewidth',2)
%         ax.XColor = 'c';
%         ax.YColor = 'c';
        set(gca,'XTick',[])
        set(gca,'YTick',[])
%         caxis([0 0.4]);
        
        annotation('rectangle',[0.65,0.415,0.05,0.07],'LineWidth',0.8)
        annotation('arrow',[0.675,0.685],[0.415,0.365],'LineWidth',1)
         
  else
        axes('Position',[0.523,0.583,0.34,0.338]);
        mesh(xx,yy,phi,'FaceColor','interp', 'EdgeColor','interp');
%         mesh(xx(10:50,10:50),yy(10:50,10:50),phi(10:50,10:50),'FaceColor','interp', 'EdgeColor','interp');
        xlim([310 390]);
        ylim([360 440]);
%         zlim([0,0.4])
        view(0,90);
        
        ax = gca;        
        set(gca,'linewidth',2)
%         ax.XColor = 'c';
%         ax.YColor = 'c';
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        box on        
        axis square;
%         caxis([0 0.4]);
        
        annotation('rectangle',[0.45,0.475,0.06,0.08],'LineWidth',0.8)
        annotation('arrow',[0.49,0.563],[0.556,0.72],'LineWidth',1)
    end

    drawnow;
    kk= kk+1;
    
%     figname = [figdir '/phi_t=' num2str(t) '.eps'];
%     print(figname,'-depsc2', '-r120')

    figname = ['../figure_MPFC_SAV/',figdir '_phi_t=' num2str(t) '.jpg'];
    print(figname,'-djpeg', '-r300')
end
