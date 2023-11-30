clear
close all

addpath(genpath('./data_sets'))
addpath(genpath('./routines'))

rng(1)
%% Set parameters
M1=10^3; M2=50; % number of data points
delta_t=0.25; % time step
ODEFUN=@(t,y) [y(2);-0.5*y(2)+y(1)-y(1).^3];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
N=2000;
PHI = @(r) exp(-r); % radial basis function used (others also work well)
x_pts = -1.2:0.01:1.2;    y_pts = -0.05:0.01:1.2;
v=(10.^(-2:0.1:1));

%% Produce the data
X = zeros(2,M1*M2);
Y=[];
for jj=1:M1
    Y0=(rand(2,1)-0.5)*4;
    [~,Y1]=ode45(ODEFUN,[0 0.000001 (1:(3+M2))*delta_t],Y0,options);
    Y1=Y1';
    X(:,(jj-1)*M2+1:jj*M2) = Y1(:,[1,3:M2+1]);
    Y(:,(jj-1)*M2+1:jj*M2) = Y1(:,3:M2+2);
end
M = M1*M2;

d=mean(vecnorm(X-mean(X,2))); % scaling for radial function

[~,C] = kmeans([X';Y'],N); % find centers

PX = zeros(M,N); PY = zeros(M,N);

for j = 1:N
    R = sqrt((X(1,:)-C(j,1)).^2+(X(2,:)-C(j,2)).^2);
    PX(:,j) = PHI(R(:)/d);
    R = sqrt((Y(1,:)-C(j,1)).^2+(Y(2,:)-C(j,2)).^2);
    PY(:,j) = PHI(R(:)/d);
end
%% EDMD

K=PX(1:M,:)\PY(1:M,:);
[V,LAM]=eig(K,'vector');
%%

I1 = find(abs(LAM-exp(1*(-0.250 + 1.392i)*delta_t))==min(abs(LAM-exp(1*(-0.250 + 1.392i)*delta_t))));
I2 = find(abs(LAM-exp(2*(-0.250 + 1.392i)*delta_t))==min(abs(LAM-exp(2*(-0.250 + 1.392i)*delta_t))));
I3 = find(abs(LAM-exp(3*(-0.250 + 1.392i)*delta_t))==min(abs(LAM-exp(3*(-0.250 + 1.392i)*delta_t))));

figure
plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'g','linewidth',2)
hold on
plot(real(LAM),imag(LAM),'.r','markersize',12)
plot(real(LAM([I1,I2,I3])),imag(LAM([I1,I2,I3])),'.b','markersize',20)
plot(real(LAM([I1,I2,I3])),-imag(LAM([I1,I2,I3])),'.b','markersize',20)
ax=gca; ax.FontSize=18; axis equal tight;   axis([x_pts(1),x_pts(end),-y_pts(end),y_pts(end)])
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
title('EDMD Eigenvalues','interpreter','latex','fontsize',18)

exportgraphics(gcf,'saved_figures/EDMD_evals_1.pdf','ContentType','vector','BackgroundColor','w')

%% Plot the eigenfunctions
[~,I]= sort(abs(LAM-1),'ascend');
V = V(:,I);

u=real([PX;PY]*V(:,2));
figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter([X(1,:),Y(1,:)],[X(2,:),Y(2,:)],5,u,'filled');
colormap(brighten(brewermap([],'RdYlBu'),-0.3))   
axis(axes1,'tight'); grid off; axis off
set(axes1,'DataAspectRatio',[1 1 1]);
pause(1)
exportgraphics(gcf,'saved_figures/EDMDefun_1.png','ContentType','image','BackgroundColor','w','Resolution',200)


%% Rerun on basin

II = find(u(1:M)>mean(u));
M = length(II);
X = X(:,II); Y = Y(:,II);
d=mean(vecnorm(X-mean(X,2))); % scaling for radial function

[~,C] = kmeans([X';Y'],N); % find centers

PX = zeros(M,N); PY = zeros(M,N);

for j = 1:N
    R = sqrt((X(1,:)-C(j,1)).^2+(X(2,:)-C(j,2)).^2);
    PX(:,j) = PHI(R(:)/d);
    R = sqrt((Y(1,:)-C(j,1)).^2+(Y(2,:)-C(j,2)).^2);
    PY(:,j) = PHI(R(:)/d);
end

K=PX(1:M,:)\PY(1:M,:);
[V,LAM]=eig(K,'vector');

I1 = find(abs(LAM-exp(1*(-0.250 + 1.392i)*delta_t))==min(abs(LAM-exp(1*(-0.250 + 1.392i)*delta_t))));
I2 = find(abs(LAM-exp(2*(-0.250 + 1.392i)*delta_t))==min(abs(LAM-exp(2*(-0.250 + 1.392i)*delta_t))));
I3 = find(abs(LAM-exp(3*(-0.250 + 1.392i)*delta_t))==min(abs(LAM-exp(3*(-0.250 + 1.392i)*delta_t))));

figure
plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'g','linewidth',2)
hold on
plot(real(LAM),imag(LAM),'.r','markersize',12)
plot(real(LAM([I1,I2,I3])),imag(LAM([I1,I2,I3])),'.b','markersize',20)
plot(real(LAM([I1,I2,I3])),-imag(LAM([I1,I2,I3])),'.b','markersize',20)
ax=gca; ax.FontSize=18; axis equal tight;   axis([x_pts(1),x_pts(end),-y_pts(end),y_pts(end)])
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
title('EDMD Eigenvalues (basin)','interpreter','latex','fontsize',18)
exportgraphics(gcf,'saved_figures/EDMD_evals_2.pdf','ContentType','vector','BackgroundColor','w')
%% Plot the eigenfunctions corresponding to lattice structure

close all

[~,I]= sort(abs(LAM-exp(1*(-0.250 + 1.392i)*delta_t)),'ascend');
LAM = LAM(I); V = V(:,I);
u=[PX;PY]*V(:,1);

figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter([X(1,:),Y(1,:)],[X(2,:),Y(2,:)],5,abs(u),'filled');
colormap(brighten(brewermap([],'RdYlBu'),-0.3))
clim([0,max(rmoutliers(abs(u)))])
axis(axes1,'tight'); grid off; axis off
set(axes1,'DataAspectRatio',[1 1 1]);
exportgraphics(gcf,'saved_figures/EDMDefun_2a.png','ContentType','image','BackgroundColor','w','Resolution',200)

figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter([X(1,:),Y(1,:)],[X(2,:),Y(2,:)],5,angle(u),'filled');
colormap(brighten(brewermap([],'RdYlBu'),-0.3))
axis(axes1,'tight'); grid off; axis off
set(axes1,'DataAspectRatio',[1 1 1]);
exportgraphics(gcf,'saved_figures/EDMDefun_2b.png','ContentType','image','BackgroundColor','w','Resolution',200)



[~,I]= sort(abs(LAM-exp(2*(-0.250 + 1.392i)*delta_t)),'ascend');
LAM = LAM(I); V = V(:,I);
u=[PX;PY]*V(:,1);


figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter([X(1,:),Y(1,:)],[X(2,:),Y(2,:)],5,abs(u),'filled');
colormap(brighten(brewermap([],'RdYlBu'),-0.3))
clim([0,max(rmoutliers(abs(u)))])
axis(axes1,'tight'); grid off; axis off
set(axes1,'DataAspectRatio',[1 1 1]);
exportgraphics(gcf,'saved_figures/EDMDefun_3a.png','ContentType','image','BackgroundColor','w','Resolution',200)

figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter([X(1,:),Y(1,:)],[X(2,:),Y(2,:)],5,angle(u),'filled');
colormap(brighten(brewermap([],'RdYlBu'),-0.3))
axis(axes1,'tight'); grid off; axis off
set(axes1,'DataAspectRatio',[1 1 1]);
exportgraphics(gcf,'saved_figures/EDMDefun_3b.png','ContentType','image','BackgroundColor','w','Resolution',200)


[~,I]= sort(abs(LAM-exp(3*(-0.250 + 1.392i)*delta_t)),'ascend');
LAM = LAM(I); V = V(:,I);
u=[PX;PY]*V(:,1);

figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter([X(1,:),Y(1,:)],[X(2,:),Y(2,:)],5,abs(u),'filled');
colormap(brighten(brewermap([],'RdYlBu'),-0.3))
clim([0,max(rmoutliers(abs(u)))])
axis(axes1,'tight'); grid off; axis off
set(axes1,'DataAspectRatio',[1 1 1]);
exportgraphics(gcf,'saved_figures/EDMDefun_4a.png','ContentType','image','BackgroundColor','w','Resolution',200)

figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter([X(1,:),Y(1,:)],[X(2,:),Y(2,:)],5,angle(u),'filled');
colormap(brighten(brewermap([],'RdYlBu'),-0.3))
axis(axes1,'tight'); grid off; axis off
set(axes1,'DataAspectRatio',[1 1 1]);
exportgraphics(gcf,'saved_figures/EDMDefun_4b.png','ContentType','image','BackgroundColor','w','Resolution',200)


[~,I]= sort(abs(LAM-exp(4*(-0.250 + 1.392i)*delta_t)),'ascend');
LAM = LAM(I); V = V(:,I);
u=[PX;PY]*V(:,1);

figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter([X(1,:),Y(1,:)],[X(2,:),Y(2,:)],5,abs(u),'filled');
colormap(brighten(brewermap([],'RdYlBu'),-0.3))
clim([0,max(rmoutliers(abs(u)))])
axis(axes1,'tight'); grid off; axis off
set(axes1,'DataAspectRatio',[1 1 1]);
exportgraphics(gcf,'saved_figures/EDMDefun_5a.png','ContentType','image','BackgroundColor','w','Resolution',200)

figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter([X(1,:),Y(1,:)],[X(2,:),Y(2,:)],5,angle(u),'filled');
colormap(brighten(brewermap([],'RdYlBu'),-0.3))
axis(axes1,'tight'); grid off; axis off
set(axes1,'DataAspectRatio',[1 1 1]);
exportgraphics(gcf,'saved_figures/EDMDefun_5b.png','ContentType','image','BackgroundColor','w','Resolution',200)















