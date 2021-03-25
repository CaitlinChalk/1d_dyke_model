clear 
close 

%directory where results are stored
res_dir = 'C:\Users\cmchalk\dyke_profiles\nonlinear_dyke_model\output\'; 

runtime = 20;
dta = 0.000944;
dt_str = strrep(num2str(dta) , '.' , '_' );

filename1 = ['Kc1' 'dt' dt_str 'ic1'];
filename2 = ['Kc1' 'dt' dt_str 'ic2'];
filename3 = ['Kc2' 'dt' dt_str 'ic2'];
%%
[z1,h1,Pe1,zf,hmax1,t] = getVariables(res_dir,filename1,dta);
[z2,h2,Pe2,zf,hmax2,t] = getVariables(res_dir,filename2,dta);
[z3,h3,Pe3,zf,hmax3,t] = getVariables(res_dir,filename3,dta);

%% plots of max h over time (convergence tests)

close all
figure(1); hold on
plot(t,hmax1,'ko-')
plot(t,hmax2,'bo-')
plot(t,hmax3,'go-')
xlabel('time','Interpreter','Latex','FontSize',12)
ylabel('$h_{max}$','Interpreter','Latex','FontSize',12)
legend(["Kc = 1, Initial condition 1","Kc = 1, Initial condition 2","Kc = 2, Initial condition 2"],'Interpreter','Latex','Location','east','FontSize',11);


%% plot h profile over time

close all
figure('Position', [1000 500 300 600]); hold on

for i = 1:length(h1)
    z = vertcat(z1{i},flipud(z1{i}),zclose);
    h = vertcat(h1{i},-flipud(h1{i}),hclose);
    plot(h,z,'k')
end

%% pressure contour plots
%see the link below for a potentially better option
%https://blogs.mathworks.com/graphics/2016/01/14/polygon-interpolation/

%times to plot
t_plot = [0,1,2,5,10,20];
[ti, idx] = min(abs(t - t_plot));

% close polygon for fill to plot properly
hclose = (-1:0.1:1)';
Peclose = zeros(length(hclose),1);
zclose = zeros(length(hclose),1);
zclose = zclose + z1{1}(1);

close all
figure('Position', [1000 500 1500 600]); hold on


for i = 1:length(t_plot)
    subplot(1,length(t_plot),i)
    z = vertcat(z1{idx(i)},flipud(z1{idx(i)}),zclose);
    h = vertcat(h1{idx(i)},-flipud(h1{idx(i)}),hclose);
    Pe = vertcat(Pe1{idx(i)},flipud(Pe1{idx(i)}),Peclose);

    p = patch(h,z,Pe);
    ylim([-100 -55]);
    xlim([-2 2]);
    c = colorbar;
    %caxis([0 1]);
    title(['t =' num2str(round(t(idx(i))))],'Interpreter','Latex','FontSize',12)
    if i == 1
        xlabel('$x/x^*$','Interpreter','Latex','FontSize',12)
        ylabel('$z/z^*$','Interpreter','Latex','FontSize',12)
    else
        set(gca,'YTick',[],'XTick',[])       
    end
end

%saveas(gcf,"../figures/pressure_contour_kc1.png")

%% make gif

close all
filename = '../figures/dyke_evolution_kc2.gif';
makegif(filename,z2,h2,Pe2,t)





