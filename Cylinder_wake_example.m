clear
close all

addpath(genpath('./data_sets'))
addpath(genpath('./routines'))

load('Cylinder_wake_data.mat')

%% Perform exact DMD
r = 47;
M = 24*5;
ind = (1:M);
X = DATA(1:160000,ind); % only look at x component of velocity
Y = DATA(1:160000,ind+1);
[~,LAM,Phi] = exactDMD(X,Y,r);
[~,I]=sort(abs(1-LAM),'ascend'); % reorder modes
Phi = Phi(:,I); LAM = LAM(I);

%% Plot the eigenvalues
figure
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
hold on
plot(real(LAM),imag(LAM),'b.','markersize',15)
axis equal
axis([-1.15,1.15,-1.15,1.15])

for j=0:5
    text(1.05*real(LAM(max(1,2*j))),1.05*imag(LAM(max(1,2*j))),sprintf('%d',j),'interpreter','latex','fontsize',13)

end
text(1.05*real(LAM(max(1,2*7)))+0.03,1.05*imag(LAM(max(1,2*7)))-0.03,'$\ddots$','interpreter','latex','fontsize',13,'rotation',-5)

title('DMD Eigenvalues','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
exportgraphics(gcf,'saved_figures/cylinder_wake_eigenvalues.pdf','ContentType','vector','BackgroundColor','none')

%% Predictive power
b = Phi\DATA(1:160000,1);
Er = zeros(1000,1);
m = mean(DATA(1:160000,1:M),2);

for j=1:1000
    Er(j)=norm(Phi*(LAM.^(j).*b)-DATA(1:160000,j+1))/norm(DATA(1:160000,j+1)-m);
end
%%
figure
loglog(Er,'b','linewidth',2)
hold on
plot([120,120],[10^(-10),max(Er)],'--k','linewidth',2)
title('Relative Prediction Error','interpreter','latex','fontsize',18)
xlabel('Number of time steps ($n$)','interpreter','latex','fontsize',18)


xx = [0.5 0.6]+0.05;    yy = [0.8 0.7];
annotation('textarrow',xx,yy,'String','Extent of snapshot data','interpreter','latex','fontsize',16,'Linewidth',1)

ax=gca; ax.FontSize=18;
exportgraphics(gcf,'saved_figures/cylinder_wake_prediction.pdf','ContentType','vector','BackgroundColor','none')


%% Plot modes
figure
tiledlayout(3,2,"TileSpacing","compact")
ct=1;
for j=[1,2,4,6,8,10]
    nexttile
    C = (reshape(Phi(1:160000,j),[800,200]));
    if mod(ct,2)==1
        C=real(C)+fliplr(real(C));
    else
        C=real(C)-fliplr(real(C));
    end
    vv=0.025:0.05:1;
    a=max(real(C(:)));
    b=min(real(C(:)));
    c=(a+b)/2;
    [~,h]=contourf((x-obst_x)/d,0.94*(y-obst_y+2.5)/d,real(C),[-vv,vv]*(b-a)+c,'edgecolor','k');
    h.LineWidth = 0.001;
    colormap(brighten(brewermap([],'RdYlBu'),0.1))
    axis equal
    hold on
    fill(obst_r*cos(0:0.01:2*pi)/d,obst_r*sin(0:0.01:2*pi)/d,'m','edgecolor','none')
    xlim([-2,10+0*max((x(:)-obst_x)/d)])
    title(sprintf('Mode %d',ct-1),'interpreter','latex','fontsize',12)
    
    pause(0.01)
    ct =ct+1;
end

exportgraphics(gcf,'saved_figures/cylinder_wake_modesx.pdf','BackgroundColor','none','Resolution',400)


