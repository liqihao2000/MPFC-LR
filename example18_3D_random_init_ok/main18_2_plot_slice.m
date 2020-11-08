clear;
clc;
close all;
dirname = 'ex18_2_MPFCdata';

datadir = [dirname,'/data'];
figdir  = [dirname,'_'];

tsave = 10;
nsave = 2000;

L = 128;
N = 128;

domain.xa  = 0;
domain.xb  = L;
domain.ya  = 0;
domain.yb  = L;
domain.za  = 0;
domain.zb  = L;

Nx = N; Ny = N; Nz=N;
hx = L/Nx;
hy = L/Ny;
hz = L/Nz;

x  = domain.xa   + hx*(0:Nx-1);
y  = domain.ya + hy*(0:Ny-1);
z  = domain.za + hy*(0:Nz-1);

[xx,yy,zz] = meshgrid(x,y,z);

figure(2)

% for t= [40:20:100 200:100:400 1000 3000]
for t= [40 80 100 300 400 3000]
    t
    filename = [datadir '/phi_t=' num2str(t)];
    ss = [filename '.txt'];
    phi = load(ss);
    n   = size(phi(:),1);
    n   = round(n^(1/3));
    phi = reshape(phi,n,n,n);
    
    xslice = []; yslice = 64; zslice = 64;
    h1 = slice(xx,yy,zz,phi, xslice,yslice,zslice);
    set(h1,'FaceColor','interp', 'EdgeColor','none')
%     colormap hsv 
    colormap jet  
    colorbar('Position',[0.845 0.35 0.03 0.52],'Fontsize',15);
    set(gca,'FontSize',18);

%     axis off;
    box on;
    axis square;
    axis equal;
    view(-51,27)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis([0 L 0 L 0 L])
    set(gca,'xtick',0:L/2:L)
    set(gca,'ytick',0:L/2:L)
    set(gca,'ztick',0:L/2:L)
    
%     figname = [figdir '/phi_t=' num2str(t) '_2.eps'];
%     print(figname,'-depsc')
    drawnow;

    figname = ['../figure_MPFC_SAV/',figdir '_phi_t=' num2str(t) '.jpg'];
    print(figname,'-djpeg', '-r300')
end
