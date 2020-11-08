clear;
clc;
close all;
dirname = 'ex19_2_1_MPFCdata';
% dirname = 'ex19_2_2_MPFCdata';
% dirname = 'ex19_2_3_MPFCdata';

datadir = [dirname,'/data'];
figdir  = [dirname,'_'];

tsave = 20;

figure(4) 
for t= 5000:tsave:5000
% for t= 2900:tsave:2900
    t
    filename = [datadir '/phi_t=' num2str(t)];
    ss = [filename '.txt'];
    phi = load(ss);
    n   = size(phi(:),1);
    n   = round(n^(1/3));
    phi = reshape(phi,n,n,n);
    
    xslice = [1 n]; yslice = [1 n]; zslice = [1 n];
    h1 = slice(phi, xslice,yslice,zslice);
    set(h1,'FaceColor','interp', 'EdgeColor','none')
%     colormap hsv
    colormap jet
    colorbar('Position',[0.815 0.30 0.03 0.55],'Fontsize',15);
    set(gca,'FontSize',18);

    axis off;
%     axis square;
    axis equal; 
    view(123,30)
    % xlabel('x')
    % ylabel('y')
    % zlabel('z')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'ztick',[])
    drawnow;
    
% %     figname = [figdir '/phi_t=' num2str(t) '_2.eps'];
% %     print(figname,'-depsc')
% 
    figname = ['../figure_MPFC_SAV/',figdir '_phi_t=' num2str(t) '_2.jpg'];
    print(figname,'-djpeg', '-r300')
end
