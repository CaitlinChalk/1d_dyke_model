function makegif(filename,z,h,Pe,t) 

% close polygon for fill to plot properly
hclose = (-1:0.1:1)';
Peclose = zeros(length(hclose),1);
zclose = zeros(length(hclose),1);
zclose = zclose + z{1}(1);

dtime = 0.01; 
for i = 1:length(h)
    %draw plot
    zp = vertcat(z{i},flipud(z{i}),zclose);
    hp = vertcat(h{i},-flipud(h{i}),hclose);
    Pep = vertcat(Pe{i},flipud(Pe{i}),Peclose);
    axis tight manual % this ensures that getframe() returns a consistent size
    f = figure('Position', [1000 500 300 600]);
    patch(hp,zp,Pep);
    ylim([-100 -55]);
    xlim([-2 2]);
    c = colorbar;
    caxis([0 1]);
    ctitle = get(c,'Title');
    set(ctitle,'String','$Pe$','Interpreter','Latex','FontSize',12);
    title(['t =' num2str(round(t(i),4))],'Interpreter','Latex','FontSize',12)
    xlabel('$x/x^*$','Interpreter','Latex','FontSize',12)
    ylabel('$z/z^*$','Interpreter','Latex','FontSize',12)
    drawnow 
      % Capture the plot as an image 
      frame = getframe(f); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,filename,'gif', 'DelayTime', dtime, 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','DelayTime', dtime, 'WriteMode','append'); 
      end 
end
close all

end