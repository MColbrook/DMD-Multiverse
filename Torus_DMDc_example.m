clear
close all

addpath(genpath('./data_sets'))
addpath(genpath('./routines'))

%% Load dataset
load('torus_data.mat')
X = X_forced(:,1:end-1);
Y = X_forced(:,2:end);
Upsilon = U(:,1:end-1);

%% Execute exactDMD and DMDc
[~,LAM,Phi] = exactDMD(X,Y,10);
w = log(LAM)/0.01; % work with logarithms of eigenvalues

[~,~,LAMc,Phic] = DMDc(X,Y,Upsilon,11,10);
wc = log(LAMc)/0.01; % work with logarithms of eigenvalues
%% Plote the (logarithm) of eigenvalues
figure
plot([imag(wTrue),-imag(wTrue)],[real(wTrue),real(wTrue)],'xb','markersize',15,'linewidth',2)
hold on
plot(imag(wc),real(wc),'.','markersize',20)
plot(imag(w),real(w),'.','markersize',20)
plot([-40,40],[0,0],'-k')
xlim([-11,11]); ylim([-0.1,0.5]);
title('Eigenvalues','interpreter','latex','fontsize',18)
xlabel('$\arg(\lambda)/\Delta t$','interpreter','latex','fontsize',18)
ylabel('$\log(|\lambda|)/\Delta t$','interpreter','latex','fontsize',18)
legend({'True','DMDc','exactDMD'},'interpreter','latex','fontsize',16,'location','north')
pbaspect([4 1.5 1])
ax=gca; ax.FontSize=18;
exportgraphics(gcf,'saved_figures/torus_eigenvalues.pdf','ContentType','vector','BackgroundColor','none')
return
%% Plot the modes
close all
clc
I = find(imag(LAM)>0);
Ic = find(imag(LAMc)>0);
for jj=1:5
    % modes defined up to phase so make phase match that of true modes
    I2 = find(abs(wTrue-wc(Ic(jj)))==min(abs(wTrue-wc(Ic(jj)))));

    uc = reshape(Phic(:,Ic(jj)),128,128);
    uTrue = reshape(PhiTrue(:,I2),128,128);
    t = uTrue./uc;
    t = mean(t(:));
    uc = t*uc;
    % subspacea(real(uc(:)),real(uTrue(:)))

    figure
    imagesc(real(uc))
    colormap(brighten(brewermap([],'RdYlBu'),0))
    axis equal off
    exportgraphics(gcf,sprintf('saved_figures/torus_DMDc_modes%d.png',jj),'BackgroundColor','none','Resolution',400)

    figure
    imagesc(real(uTrue))
    colormap(brighten(brewermap([],'RdYlBu'),0))
    axis equal off
    exportgraphics(gcf,sprintf('saved_figures/torus_true_modes%d.png',jj),'BackgroundColor','none','Resolution',400)
    
    find(abs(w-wc(Ic(jj)))==min(abs(w-wc(Ic(jj)))))
    u = reshape(Phi(:,abs(w-wc(Ic(jj)))==min(abs(w-wc(Ic(jj))))),128,128);
    
    t = uTrue./u;
    t = mean(t(:));
    u = t*u;
    % subspacea(real(u(:)),real(uTrue(:)))
    figure
    imagesc(real(u))
    colormap(brighten(brewermap([],'RdYlBu'),0))
    axis equal off
    exportgraphics(gcf,sprintf('saved_figures/torus_DMD_modes%d.png',jj),'BackgroundColor','none','Resolution',400)

    u = reshape(Phi(:,I(jj)),128,128);
    imagesc(real(u))
    colormap(brighten(brewermap([],'RdYlBu'),0))
    axis equal off
    exportgraphics(gcf,sprintf('saved_figures/torus_DMD_modes_unordered%d.png',jj),'BackgroundColor','none','Resolution',400)
    close all

end

