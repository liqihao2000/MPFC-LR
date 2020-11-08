clear;
clc;
close all;
dirname = 'ex19_2_1_MPFCdata';
% dirname = 'ex19_2_2_MPFCdata';
% dirname = 'ex19_2_3_MPFCdata';

datadir = [dirname,'/data'];
figdir  = [dirname,'_'];

tsave = 20;

figure(3)
for t= 5000:tsave:5000
% for t= 2980:tsave:3000
    t
    filename = [datadir '/phi_t=' num2str(t)];
    ss = [filename '.txt'];
    phi = load(ss);
    n   = size(phi(:),1);
    n   = round(n^(1/3));
    phi = reshape(phi,n,n,n);
    
    p1 = patch(isosurface(phi,0));
    isonormals(phi,p1);
    set(p1,'FaceColor','y', 'EdgeColor','none');daspect([1 1 1]); 
    camlight;  
    lighting phong;
    box on; axis image;
    %    axis([min(xx(:)) max(xx(:)) min(yy(:)) max(yy(:)) min(zz(:)) max(zz(:)) ]);
    % view(-161,16);
    view(123,30)
    % xlabel('x')
    % ylabel('y')
    % zlabel('z')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'ztick',[])
    drawnow;
    
% %     figname = [figdir '/phi_t=' num2str(t) '.eps'];
% %     print(figname,'-depsc')

    figname = ['../figure_MPFC_SAV/',figdir '_phi_t=' num2str(t) '.jpg'];
    print(figname,'-djpeg', '-r300')
end
