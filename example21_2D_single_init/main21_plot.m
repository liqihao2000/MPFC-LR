clear;
clc;
close all;
dirname = 'ex21_1_MPFCdata';

datadir = [dirname,'/data'];
figdir  = [dirname,'_'];

tsave = 100;
nsave = 1000;

Nx = 512; Ny = Nx;
domain.left   = 0; 
domain.right  = 340;
domain.bottom = 0;
domain.top    = 340;
Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;
hx = Lx/Nx;
hy = Ly/Ny;

x  = domain.left   + hx*(0:Nx-1);
y  = domain.bottom + hy*(0:Ny-1);

[xx,yy] = meshgrid(x,y);

kk=1;

% for t= 0:tsave:nsave
% for t= [ 0]
% for t= [ 0 100 300 500 700 900 1100 3000]
for t= [ 0 500 1100 3000]
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


        axes('Position',[0.603,0.673,0.25,0.25]);
        mesh(xx,yy,phi,'FaceColor','interp', 'EdgeColor','interp');
%         mesh(xx(10:50,10:50),yy(10:50,10:50),phi(10:50,10:50),'FaceColor','interp', 'EdgeColor','interp');
        xlim([170-30 170+30]);
        ylim([170-30 170+30]);
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
        
        annotation('rectangle',[0.47,0.463,0.09,0.115],'LineWidth',0.8)
        annotation('arrow',[0.52,0.632],[0.58,0.78],'LineWidth',1)

    
    drawnow;
    kk = kk +1;
    
%     figname = [figdir '/phi_t=' num2str(t) '.eps'];
%     print(figname,'-depsc2', '-r120')

    figname = ['../figure_MPFC_SAV/',figdir '_phi_t=' num2str(t) '.jpg'];
    print(figname,'-djpeg', '-r300')
end
