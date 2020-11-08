function showsolution_3D(nfigure,xx,yy,zz,phi0,t,dir_fig,format_fig)
% figure(nfigure);
xslice = []; yslice = 25; zslice = 25;
h1 = slice(xx,yy,zz,phi0, xslice,yslice,zslice);
set(h1,'FaceColor','interp', 'EdgeColor','none')
% camproj perspective;
box on;
% set(h1,'looseInset',[0 0 0 0])
colormap jet;
colorbar;
axis square;
% axis tight;
% % axis off;
axis equal;  
view(-55,23)
% caxis( [-max(abs(phi0(:))) max(abs(phi0(:))) ]);
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
