clear
close all

addpath(genpath('./data_sets'))
addpath(genpath('./routines'))

load('rossler_data.mat')
% This data was produced using the code of Claire Valva: https://github.com/clairevalva/resolvent_minimal
% which also uses the code of Dimitrios Giannakis: https://github.com/dg227/NLSA

%% Plot the approximate eigenfunctions and various trajectories
T = round(150/dt);

for IND = 1:3
    example_eigenfs(:,IND)=example_eigenfs(:,IND)/mean(abs(example_eigenfs(1:T,IND)));

    figure
    u = imag(example_eigenfs(:,IND));
    scatter3(orig_data(1,:),orig_data(2,:),orig_data(3,:),7,u,'filled');
    colormap(brighten(inferno,0.3))
    colormap(brighten(brewermap([],'RdYlBu'),-0.4))
    clim([-3*mean(abs(u)),3*mean(abs(u))])
    axis tight; grid off; axis off
    axis equal
    set(gca,'DataAspectRatio',[1 1 1]);
    xlabel('$X$','interpreter','latex','fontsize',14)
    ylabel('$Y$','interpreter','latex','fontsize',14)
    zlabel('$Z$','interpreter','latex','fontsize',14,'rotation',0)
    view(gca,[3.7779   17.9883]);
    if IND ==1
        title({'$\sigma\approx 1.0261$'},'interpreter','latex','fontsize',18)
    elseif IND==2
        title({'$\sigma\approx 2.0516$'},'interpreter','latex','fontsize',18)
    else
        title({'$\sigma\approx 0.3785$'},'interpreter','latex','fontsize',18)
    end
    exportgraphics(gcf,sprintf('saved_figures/rossler_efun_%d.png',IND),'ContentType','image','BackgroundColor','w','Resolution',200)
    
    figure
    plot((0:length(u)-1)*dt,real(example_eigenfs(:,IND)),'b','linewidth',2)
    xlim([0,round(T*dt)])
    xlabel('Time','interpreter','latex','fontsize',18)
    title('$\mathrm{Re}(\phi)$','interpreter','latex','fontsize',18)
    ax=gca; ax.FontSize=18;
    exportgraphics(gcf,sprintf('saved_figures/rossler_wave_%d.png',IND),'ContentType','image','BackgroundColor','w','Resolution',200)
    
    
    figure
    plot(real(example_eigenfs(1:T,IND)),imag(example_eigenfs(1:T,IND)),'b','linewidth',2)
    hold on
    plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'g','linewidth',3)
    axis equal
    axis([-1,1,-1,1]*1.6)
    % title('DMD Eigenvalues','interpreter','latex','fontsize',18)
    xlabel('$\mathrm{Re}(\phi)$','interpreter','latex','fontsize',18)
    ylabel('$\mathrm{Im}(\phi)$','interpreter','latex','fontsize',18)
    ax=gca; ax.FontSize=18;
    exportgraphics(gcf,sprintf('saved_figures/rossler_spiral_%d.png',IND),'ContentType','image','BackgroundColor','w','Resolution',200)
    close all
end



%% Residual plot
clear
load('rossler_data.mat')
T = round(1000/dt);
res = zeros(T,3);

for IND = 1:3
    example_eigenfs(:,IND)=example_eigenfs(:,IND)/norm(example_eigenfs(:,IND));
    lam=example_freqs(IND);
    for del = 1:T
        u1 = example_eigenfs(1:end-del,IND);
        u2 = example_eigenfs((1+del):end,IND);
        res(del,IND) = norm(u2-exp(1i*lam*dt*del)*u1)/norm(u2);
    end
end
%%
figure
loglog((1:T)*dt,res,'linewidth',2)
hold on
plot([1,1]/0.071,[0.001,10],'--k','linewidth',2)
xlabel('Time (s)','interpreter','latex','fontsize',18)
title('$\|\phi(t)-\exp(i\sigma t)\phi(0)\|/\|\phi\|$','interpreter','latex','fontsize',18)
legend({'$\sigma\approx 1.0261$','$\sigma\approx 2.0516$','$\sigma\approx 0.3785$','Lyapunov Time'},'interpreter','latex','fontsize',16,'location','northwest')
ax=gca; ax.FontSize=18;
xlim([0.01,1000])
ylim([min(res(:)/1.1),2])
exportgraphics(gcf,'saved_figures/rossler_residuals.pdf','ContentType','vector','BackgroundColor','none')


return


