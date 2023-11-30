clear
% close all
rng(1); % set random seed to get identical snapshots each time
addpath(genpath('./routines'))

% NB general code for mpEDMD can be found at: https://github.com/MColbrook/Measure-preserving-Extended-Dynamic-Mode-Decomposition

%% Set parameters
options = odeset('RelTol',1e-14,'AbsTol',1e-14); % for the numerical solver
SIGMA=10;   BETA=8/3;   RHO=28;
ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));y(1).*(RHO-y(3))-y(2);y(1).*y(2)-BETA*y(3)];

N=500;                              % number of delay embeddings
g = @(x,y,z) tanh((x.*y-3*z)/5);    % observable
M=10^4;                             % number of data points
dt=0.1;                             % time step for trajectory sampling
  
%% Produce the data - slowest part of this script!
Y0=(rand(3,1)-0.5)*4;
[~,Y0]=ode45(ODEFUN,[0.000001 1, 100],Y0,options); Y0 = Y0(end,:)'; % sample after when on the attractor
[~,DATA]=ode45(ODEFUN,[0.000001 (1:((M+(N+1))))*dt],Y0,options);

%% Use delay embedding
PX1=zeros(M,N); PX1(:,1)=DATA(1:M,1); PX2=zeros(M,N); PX2(:,1)=DATA(1:M,2); PX3=zeros(M,N); PX3(:,1)=DATA(1:M,3);
PY1=zeros(M,N); PY1(:,1)=DATA((1:M)+1,1); PY2=zeros(M,N); PY2(:,1)=DATA((1:M)+1,2); PY3=zeros(M,N); PY3(:,1)=DATA((1:M)+1,3);
for j=2:N
    PX1(:,j)=DATA((1:M)+(j-1),1); PX2(:,j)=DATA((1:M)+(j-1),2); PX3(:,j)=DATA((1:M)+(j-1),3);
    PY1(:,j)=DATA((1:M)+1+(j-1),1); PY2(:,j)=DATA((1:M)+1+(j-1),2); PY3(:,j)=DATA((1:M)+1+(j-1),3);
end

%% Run mpEDMD
PX = g(PX1,PX2,PX3);
PY = g(PY1,PY2,PY3);
[~,mpV,mpD] = mpEDMDqr(PX,PY,1/M);

% %% Plot the singular values
% S = svd(PX);
% figure
% semilogy(S,'b.','markersize',15)
% title('Singular Values','interpreter','latex','fontsize',18)
% ax=gca; ax.FontSize=18;
% exportgraphics(gcf,'saved_figures/Lorenz_SVD.pdf','ContentType','vector','BackgroundColor','none')

%% Compute spectral measure
c = zeros(N,1); c(1) = 1; % coefficients of g in Krylov basis
piE = diag(mpD); TH=angle(piE*exp(1i*eps)); % eps here is to take into account MATLAB's convention for angle
G = (PX'*PX)/M;
MU=abs(mpV'*G*c).^2;

%% Cdf plots
[~,Ib] = sort(TH(:),'ascend');
THp=TH(Ib); THp=[THp(:)-10^(-14),THp(:)]'; THp=THp(:);
cdf=0*THp;
cc=0;
for j=1:length(TH)
    cdf(2*j-1)=cc;    cc=cc+MU(Ib(j));    cdf(2*j)=cc;
end

THp = [-pi;THp(:);pi]; cdf = [0;cdf(:);sum(MU)]; % for visualisation


figure
plot(THp,cdf/sum(MU),'b','linewidth',2)
ylim([0,1]); xlim([-pi,pi]);

ax=gca; ax.FontSize=18;
title(sprintf('$F_{\\mu_g^{(%d)}}(\\theta)$',N),'interpreter','latex','fontsize',30)
xlabel('$\theta$','interpreter','latex','fontsize',30)
exportgraphics(gcf,sprintf('saved_figures/lorenz_cdf_%d.pdf',N),'ContentType','vector','BackgroundColor','none')
