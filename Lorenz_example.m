clear
close all

addpath(genpath('./data_sets'))
addpath(genpath('./routines'))

%% Set parameters for numerical simulation
options = odeset('RelTol',1e-14,'AbsTol',1e-14);
SIGMA=10;   BETA=8/3;   RHO=28;
ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));y(1).*(RHO-y(3))-y(2);y(1).*y(2)-BETA*y(3)];

%% Set parameters

for N=[10:10:50,100] % number of delay embeddings
    rng(1)
    M=5*10^5;     % number of data points
    dt=0.001;    % time step for trajectory sampling
    dt2=0.2;    % time step for delays
    h=dt2/dt;
    
    %% Produce the data
    Y0=(rand(3,1)-0.5)*4;
    [~,Y0]=ode45(ODEFUN,[0.000001 1, 100],Y0,options); Y0 = Y0(end,:)'; % sample after when on the attractor
    [~,DATA]=ode45(ODEFUN,[0.000001 (1:((M+h*(N+1))))*dt],Y0,options);
    
    %% Use delay embedding
    PX=zeros(M,N); PX(:,1)=DATA(1:M,1);
    PX2=zeros(M,N); PX2(:,1)=DATA(1:M,2);
    PX3=zeros(M,N); PX3(:,1)=DATA(1:M,3);
    PY=zeros(M,N); PY(:,1)=DATA((1:M)+1,1);
    PY2=zeros(M,N); PY2(:,1)=DATA((1:M)+1,2);
    PY3=zeros(M,N); PY3(:,1)=DATA((1:M)+1,3);
    
    for j=2:N
        PX(:,j)=DATA((1:M)+h*(j-1),1);
        PX2(:,j)=DATA((1:M)+h*(j-1),2);
        PX3(:,j)=DATA((1:M)+h*(j-1),3);
        PY(:,j)=DATA((1:M)+1+h*(j-1),1);
        PY2(:,j)=DATA((1:M)+1+h*(j-1),2);
        PY3(:,j)=DATA((1:M)+1+h*(j-1),3);
    end

    K = [PX,PX2,PX3]\[PY,PY2,PY3];
    [V,LAM] = eig(K,'vector');
    [~,I]=sort(abs(angle(LAM)+0.006),'ascend');   
    V =V(:,I); LAM =LAM(I);

    w = log(LAM)/dt;

    figure
    plot([-40,40],[0,0],'-k')
    hold on
    plot(imag(w),real(w),'b.','markersize',15)
    title('DMD Eigenvalues','interpreter','latex','fontsize',18)
    xlabel('$\arg(\lambda)/\Delta t$','interpreter','latex','fontsize',18)
    ylabel('$\log(|\lambda|)/\Delta t$','interpreter','latex','fontsize',18)
    ax=gca; ax.FontSize=18;
    ylim([-0.2,0.2]); xlim([-25,25])
    pbaspect([4 1.5 1])
    text(-20,0.1,sprintf('$N=%d$',N),'interpreter','latex','fontsize',18)
    exportgraphics(gcf,sprintf('saved_figures/lorenz_eigenvalues%d.pdf',N),'ContentType','vector','BackgroundColor','none')    
    
    for j=1
        figure
        u = [PX,PX2,PX3]*V(:,1);
        u = -real(u(1:min(10^6,M),:));
        u = u/norm(u);
        [~,I]=sort(u,'ascend');
        scatter3(PX(I,1),PX2(I,1),PX3(I,1),6,u(I),'filled');
        colormap(brighten(brewermap([],'RdYlBu'),-0.3))
        clim([-3*mean(abs(u)),3*mean(abs(u))])
        view(gca,[13.1786087602293 -1.28469255513244]);
        axis tight; grid off; axis off
        axis equal
        set(gca,'DataAspectRatio',[1 1 1]);
        xlabel('$X$','interpreter','latex','fontsize',14)
        ylabel('$Y$','interpreter','latex','fontsize',14)
        zlabel('$Z$','interpreter','latex','fontsize',14,'rotation',0)
        text(-5,0,43,sprintf('$N=%d$',N),'interpreter','latex','fontsize',18)

        exportgraphics(gcf,sprintf('saved_figures/lorenz_efun_angle_%d_N%d.png',j,N),'ContentType','image','BackgroundColor','w','Resolution',200)
        pause(1)
    end
end


