function showsolution_3D_2(nfigure,xx,yy,zz,phi0,t,dir_fig,format_fig)
% figure(nfigure);
p1 = patch(isosurface(xx,yy,zz,phi0,0));
set(p1,'FaceColor','g', 'EdgeColor','none');daspect([1 1 1]);
% camlight;
% lighting phong; 
box on; axis image;
axis([min(xx(:)) max(xx(:)) min(yy(:)) max(yy(:)) min(zz(:)) max(zz(:)) ]);
% view(-161,16);
view(60,10)
% xlabel('x')
% ylabel('y')
% zlabel('z')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])

% % camproj perspective;
% set(p1,'looseInset',[0 0 0 0])
% colormap jet;
% colorbar;
% axis square;
% % axis tight;
% axis off;
% axis equal;  
% % caxis( [-max(abs(phi0(:))) max(abs(phi0(:))) ]);

if exist('t','var')
    title(['$t=$',num2str(t)],'Interpreter','LaTex');
    colorbar;
end
if exist('dir_fig','var')
    if ~exist(dir_fig,'dir')
       mkdir(dir_fig); 
    end
    datafile = [dir_fig '/phi_t=' num2str(t) '.mat'];
%     save(datafile,'phi0','t');
    
    title([]);
    colorbar off;
    if exist('format_fig','var')
        figfile =  [dir_fig '/phi_t=' num2str(t) format_fig];
    else
        figfile =  [dir_fig '/phi_t=' num2str(t) '.jpg'];        
    end
    saveas(1,figfile);
end
drawnow;

end
