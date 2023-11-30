clear
close all

addpath(genpath('./saved_data_from_runs'))

% NB general code for ResDMD can be found at: https://github.com/MColbrook/Residual-Dynamic-Mode-Decomposition

%% UNCOMMENT TO RUN
% rng(1)
% %% Set parameters
% M1=10^3; % number of data points
% M2=50;
% delta_t=0.25; % time step
% ODEFUN=@(t,y) [y(2);y(1)-y(1).^3];
% options = odeset('RelTol',1e-12,'AbsTol',1e-12);
% N=100;
% PHI = @(r) exp(-r); % radial basis function
% 
% x_pts = -1.2:0.04:1.2;    y_pts = -0.05:0.04:1.2;
% v=(10.^(-2:0.1:1));
% 
% %% Produce the data
% X=[];
% Y=[];
% for jj=1:M1
%     Y0=(rand(2,1)-0.5)*4;
%     [~,Y1]=ode45(ODEFUN,[0 0.000001 (1:(3+M2))*delta_t],Y0,options);
%     Y1=Y1';
%     X = [X,Y1(:,[1,3:M2+1])];
%     Y = [Y,Y1(:,3:M2+2)];
% end
% M = M1*M2;
% 
% d=mean(vecnorm(X-mean(X')')); % scaling for radial function
% 
% [~,C] = kmeans([X';Y'],N); % find centers
% 
% PX = zeros(M,N); PY = zeros(M,N);
% 
% for j = 1:N
%     R = sqrt((X(1,:)-C(j,1)).^2+(X(2,:)-C(j,2)).^2);
%     PX(:,j) = PHI(R(:)/d);
%     R = sqrt((Y(1,:)-C(j,1)).^2+(Y(2,:)-C(j,2)).^2);
%     PY(:,j) = PHI(R(:)/d);
% end
% %% EDMD
% K=PX(1:M,:)\PY(1:M,:);
% [V,LAM]=eig(K,'vector');
% 
% %% ResDMD for pseudospectrum, code available here: https://github.com/MColbrook/Residual-Dynamic-Mode-Decomposition
% 
% z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);
% 
% RES = KoopPseudoSpecQR(PX,PY,1/M,z_pts,'Parallel','off');
% RES = reshape(RES,length(y_pts),length(x_pts));

%%
close all
clear;load('duffing3.mat')

figure
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(max(min(v),real(RES))),log10(v));
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(max(min(v),real(RES))),log10(v));
cbh=colorbar;
cbh.Ticks=log10([0.005,0.01,0.1,1]);
cbh.TickLabels=[0,0.01,0.1,1];
clim([log10(min(v)),0]);
reset(gcf)
set(gca,'YDir','normal')
colormap gray
ax=gca; ax.FontSize=18; axis equal tight;   axis([x_pts(1),x_pts(end),-y_pts(end),y_pts(end)])
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
title(sprintf('$\\mathrm{Sp}_\\epsilon(\\mathcal{K})$, $N=%d$',N),'interpreter','latex','fontsize',18)

exportgraphics(gcf,sprintf('saved_figures/duffing_pseudospec_N%d.png',N),'ContentType','image','BackgroundColor','w','Resolution',200)

figure
plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'g','linewidth',2)
hold on
plot(real(LAM),imag(LAM),'.r','markersize',12)
ax=gca; ax.FontSize=18; axis equal tight;   axis([x_pts(1),x_pts(end),-y_pts(end),y_pts(end)])
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
title(sprintf('DMD Eigenvalues, $N=%d$',N),'interpreter','latex','fontsize',18)

exportgraphics(gcf,sprintf('saved_figures/duffing_evals_N%d.pdf',N),'ContentType','vector','BackgroundColor','w')

res = (vecnorm(PY*V-PX*V*diag(LAM))./vecnorm(PX*V))'; % ResDMD for residuals

figure
histogram(res,10,'FaceColor',[0.8500 0.3250 0.0980])
xlabel('$\|(\mathcal{K}-\lambda_j I)g_j\|$','interpreter','latex','fontsize',18)
title(sprintf('Histogram of Errors, $N=%d$',N),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
exportgraphics(gcf,sprintf('saved_figures/duffing_hist_N%d.pdf',N),'ContentType','vector','BackgroundColor','w')

%% Plot the eigenfunctions
[~,I]= sort(abs(1-LAM),'ascend');
LAM = LAM(I); V = V(:,I);
close all

xvec = (-2:0.25:0) - 0.01;
xvec = [xvec,-xvec];

for j=1:14

    u=real([PX;PY]*V(:,j));
    figure1=figure;
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    scatter([X(1,:),Y(1,:)],[X(2,:),Y(2,:)],5,u,'filled');
    colormap(brighten(brewermap([],'RdYlBu'),-0.3))
    clim([min(u),max(u)])
    
    for jj=1:length(xvec)
        ODEFUN=@(t,y) [y(2);y(1)-y(1).^3];
        [~,Y1]=ode45(ODEFUN,0:0.01:20,[xvec(jj),0],options);
        plot(Y1(:,1),Y1(:,2),'k','linewidth',1);
    end
    
    axis(axes1,'tight'); grid off; axis off%hold(axes1,'off');
    set(axes1,'DataAspectRatio',[1 1 1]);
    pause(1)
    exportgraphics(gcf,sprintf('saved_figures/duffing_efun_%d.png',j),'ContentType','image','BackgroundColor','w','Resolution',200)
    close all
end












