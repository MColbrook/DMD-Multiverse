clear
close all

addpath(genpath('./data_sets'))
addpath(genpath('./routines'))
addpath(genpath('./saved_data_from_runs'))

options = odeset('RelTol',1e-14,'AbsTol',1e-14);

%% Set parameters

N = 10;
Mvec = round(10.^(3:0.2:6));

% Er = zeros(length(Mvec),3,50);
% rng(1)
% 
% for S = 1:size(Er,3)
%     %% Produce the data
%     dt=0.001;    % time step for trajectory sampling
%     dt2=0.2;    % time step for delays
%     h=dt2/dt;
% 
%     SIGMA=10;   BETA=8/3;   RHO=28;
%     ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));y(1).*(RHO-y(3))-y(2);y(1).*y(2)-BETA*y(3)];
%     Y0=(rand(3,1)-0.5)*4;
%     [~,Y0]=ode45(ODEFUN,[0.000001 1, 100],Y0,options); Y0 = Y0(end,:)'; % sample after when on the attractor
%     [~,DATA]=ode45(ODEFUN,[0.000001 (1:((max(Mvec)+h*(N+1))))*dt],Y0,options);
%     ctt = 1;
%     for M = Mvec  % number of delay embeddings
% 
% 
% 
% 
% 
% 
%         %% Use delay embedding
%         PX=zeros(M,N); PX(:,1)=DATA(1:M,1);
%         PX2=zeros(M,N); PX2(:,1)=DATA(1:M,2);
%         PX3=zeros(M,N); PX3(:,1)=DATA(1:M,3);
%         PY=zeros(M,N); PY(:,1)=DATA((1:M)+1,1);
%         PY2=zeros(M,N); PY2(:,1)=DATA((1:M)+1,2);
%         PY3=zeros(M,N); PY3(:,1)=DATA((1:M)+1,3);
% 
%         for j=2:N
%             PX(:,j)=DATA((1:M)+h*(j-1),1);
%             PX2(:,j)=DATA((1:M)+h*(j-1),2);
%             PX3(:,j)=DATA((1:M)+h*(j-1),3);
%             PY(:,j)=DATA((1:M)+1+h*(j-1),1);
%             PY2(:,j)=DATA((1:M)+1+h*(j-1),2);
%             PY3(:,j)=DATA((1:M)+1+h*(j-1),3);
%         end
%         %%
%         G = (PX'*PX)/M;
%         A = (PX'*PY)/M;
%         L = (PY'*PY)/M;
%         G=([PX,PX2,PX3])'*([PX,PX2,PX3])/M;
%         A=([PX,PX2,PX3])'*([PY,PY2,PY3])/M;
%         K = [PX,PX2,PX3]\[PY,PY2,PY3];
%         % 
%         % [~,V,LAM] = mpEDMD(G,A); LAM=diag(LAM);LAM=LAM(:);
%         % [~,I]=sort(abs((LAM)),'descend');
% 
%         % 
% 
%         % 
%         [~,LAM] = eig(K,'vector');
%         w = log(LAM)/dt;
%         Er(ctt,1,S)=mean(abs(real(w)));
% 
% 
% 
%         % [~,LAM] = fbDMD(([PX,PX2,PX3])',([PY,PY2,PY3])',size(PX,2)*3);
%         % w = log(LAM)/dt;
%         % Er(ctt,2)=mean(abs(real(w)));
% 
%         [~,LAM] = tlsDMD(([PX,PX2,PX3])',([PY,PY2,PY3])',size(PX,2)*3);
%         w = log(LAM)/dt;
%         Er(ctt,2,S)=mean(abs(real(w)));
% 
% 
%         [~,LAM] = optdmd(([PX,PX2,PX3])',(1:M)*dt,3*N,2,varpro_opts('maxiter',100),1i*imag(log(LAM)/dt));%varpro_opts('ifprint',0,'maxiter',10,'tol',10^(-9))
%         w = LAM;
%         Er(ctt,3,S)=mean(abs(real(w)));
%     close all
%     loglog(Mvec,mean(squeeze(Er(:,1,:)),2),'k');
%     hold on
%     loglog(Mvec,mean(squeeze(Er(:,2,:)),2),'g');
%     loglog(Mvec,mean(squeeze(Er(:,3,:)),2),'b');
%     pause(0.1)
% 
%         ctt = ctt+1;
%     end
% 
% 
% 
% 
% end

%% Plot the figure
load('Lorenz_Denoise.mat')
figure
loglog(Mvec,Er(:,1),'.-k','linewidth',1,'markersize',12)
hold on
loglog(Mvec,Er(:,2),'.-r','linewidth',1,'markersize',12)
loglog(Mvec,Er(:,3),'.-b','linewidth',1,'markersize',12)
legend({'exactDMD','tlsDMD','optDMD'},'interpreter','latex','fontsize',12,'location','southwest')
title('Mean $|\log(|\lambda|)/\Delta t|$','interpreter','latex','fontsize',18)
xlabel('$M$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
% ylim([10^(-15),1])
exportgraphics(gcf,'lorenz_denoised.pdf','ContentType','vector','BackgroundColor','none')






