clear
close all

addpath(genpath('./data_sets'))
addpath(genpath('./routines'))

load('Cylinder_wake_data.mat')
%% Set parameters
M = 24*5;
pvec = 1:40; % sparsity parameter for convergence plot
X = DATA(1:800*200,1:M);
Y = DATA(1:800*200,2:(M+1));

%% Run standard DMD for reference
[K,LAM,Phi] = exactDMD(X,Y,47);
[~,I]  = sort(abs(1-LAM),'ascend');
Phi = Phi(:,I); LAM = LAM(I);

% %% Run exactDMD, cDMD and rDMD - uncomment to run
% 
% rng(1);
% Er_lam = zeros(6,length(pvec),25,3);
% Er_vec = zeros(6,length(pvec),25,3);
% tm = zeros(length(pvec),1,3);
% % profile on
% for s = 1:size(Er_lam,3)
% 
%     ct = 1; % counting
%     for p = pvec
% 
%         tic
%         [~,LAM1,Phi1] = exactDMD(X,Y,p);
%         tm(ct,s,1) = toc;
%         [~,I]  = sort(abs(1-LAM1),'ascend');
%         Phi1 = Phi1(:,I); LAM1 = LAM1(I);
% 
%         Er_lam(1,ct,s,1) = min(abs(LAM(1)-LAM1));
%         Er_vec(1,ct,s,1) = subspacea(Phi(:,1),Phi1(:,1));
%         for jj=1:5
%             if ct>2*jj
%                 Er_lam(jj+1,ct,s,1) = min(abs(LAM(2*jj)-LAM1));
%                 Er_vec(jj+1,ct,s,1) = 5;
%                 for kk = 1:min(2*jj+2,p)
%                     Er_vec(jj+1,ct,s,1) = min(Er_vec(jj+1,ct,s,1),subspacea(Phi(:,2*jj),Phi1(:,kk)));
%                 end
%             end
%         end
% 
%         tic
%         [C,theta] = rand_meas(p,800,200);
%         Xc = C*X;
%         Yc = C*Y;
%         [W1,LAM1,~,~,V1,Sinv1] = standard_DMD(Xc,Yc,p);
%         Phi1 = Y*V1*Sinv1*W1;
%         tm(ct,s,2) = toc;
%         [~,I]  = sort(abs(1-LAM1),'ascend');
%         Phi1 = Phi1(:,I); LAM1 = LAM1(I);
% 
%         Er_lam(1,ct,s,2) = min(abs(LAM(1)-LAM1));
%         Er_vec(1,ct,s,2) = subspacea(Phi(:,1),Phi1(:,1));
%         for jj=1:5
%             if ct>2*jj
%                 Er_lam(jj+1,ct,s,2) = min(abs(LAM(2*jj)-LAM1));
%                 Er_vec(jj+1,ct,s,2) = 5;
%                 for kk = 1:min(2*jj+2,p)
%                     Er_vec(jj+1,ct,s,2) = min(Er_vec(jj+1,ct,s,2),subspacea(Phi(:,2*jj),Phi1(:,kk)));
%                 end
%             end
%         end
% 
%         [~,LAM1,tm(ct,s,3),Phi1] = rDMD(X,Y,p,10);
%         [~,I]  = sort(abs(1-LAM1),'ascend');
%         Phi1 = Phi1(:,I); LAM1 = LAM1(I);
% 
%         Er_lam(1,ct,s,3) = min(abs(LAM(1)-LAM1));
%         Er_vec(1,ct,s,3) = subspacea(Phi(:,1),Phi1(:,1));
%         for jj=1:5
%             if ct>2*jj
%                 Er_lam(jj+1,ct,s,3) = min(abs(LAM(2*jj)-LAM1));
%                 Er_vec(jj+1,ct,s,3) = 5;
%                 for kk = 1:min(2*jj+2,p)
%                     Er_vec(jj+1,ct,s,3) = min(Er_vec(jj+1,ct,s,3),subspacea(Phi(:,2*jj),Phi1(:,kk)));
%                 end
%             end
%         end
%     end
% 
% end


%%
load('Cylinder_Timings.mat')
close all
figure;
semilogy(pvec,mean(mean(Er_lam(:,:,:,2),3),1),'linewidth',2)
hold on
semilogy(pvec,mean(mean(Er_vec(:,:,:,2),3),1),'linewidth',2)
legend({'Eigenvalue Error','DMD Mode Error'},'interpreter','latex','fontsize',16,'location','southwest')
title('cDMD','interpreter','latex','fontsize',18)
xlabel('$p$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
ylim([10^(-15),1])
exportgraphics(gcf,'saved_figures/cylinder_compressed_error1.pdf','ContentType','vector','BackgroundColor','none')

figure;
semilogy(pvec,mean(mean(Er_lam(:,:,:,3),3),1),'linewidth',2)
hold on
semilogy(pvec,mean(mean(Er_vec(:,:,:,3),3),1),'linewidth',2)
legend({'Eigenvalue Error','DMD Mode Error'},'interpreter','latex','fontsize',16,'location','southwest')
title('rDMD','interpreter','latex','fontsize',18)
xlabel('$r$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
ylim([10^(-15),1])
exportgraphics(gcf,'saved_figures/cylinder_compressed_error2.pdf','ContentType','vector','BackgroundColor','none')

figure
loglog(mean(tm(:,:,1),2),mean(mean(Er_lam(:,:,:,1),3),1),'.','markersize',14)
hold on
loglog(mean(tm(:,:,2),2),mean(mean(Er_lam(:,:,:,2),3),1),'.','markersize',14)
loglog(mean(tm(:,:,3),2),mean(mean(Er_lam(:,:,:,3),3),1),'.','markersize',14)
legend({'exactDMD','cDMD','rDMD'},'interpreter','latex','fontsize',16,'location','southwest')
xlabel('Time (s)','interpreter','latex','fontsize',18)
ylabel('Error','interpreter','latex','fontsize',18)
grid on
ax=gca; ax.FontSize=18;
exportgraphics(gcf,'saved_figures/cylinder_compressed_error_times.pdf','ContentType','vector','BackgroundColor','none')

%% Compressed sensing DMD

rng(1);
ct = 1; % counting
for p = 800
    [C,theta] = rand_meas(p,800,200);
    Xc = C*X;
    Yc = C*Y;
    [W1,LAM1,~,~,V1,Sinv1] = standard_DMD(Xc,Yc,47);
    [~,I]  = sort(abs(1-LAM1),'ascend');
    W1 = W1(:,I); LAM1 = LAM1(I);
    Phi2 = zeros(size(Phi,1),6);
    Phic = Yc*V1*Sinv1*W1;
    for jj=1:6
        if 2*jj-1<=p
            z=CoSaMP( theta, Phic(:,2*jj-1),50, []);
            z2 = reshape(ifft2(ifftshift(reshape(z,800,200))),800*200,1)*sqrt(800*200);
            Phi2(:,jj) = z2;

        end
    end
    ct=ct+1
end

%% Plot reconstructed modes

for jj=1:6
    C = real(reshape(Phi2(:,jj),800,200));
    if jj==2 % just to make color schemes match
        C=-C;
    elseif jj==3
        C=-C;
    elseif jj==6
        C=-C;
    end
    vv=0.025:0.05:1;
    a=min(real(C(:)));
    b=max(real(C(:)));
    c=(a+b)/2;
    figure
    [~,h]=contourf((x-obst_x)/d,0.94*(y-obst_y)/d,real(C),[-vv,vv]*(b-a)+c,'edgecolor','none');
    h.LineWidth = 0.001;
    colormap(brighten(brewermap([],'RdYlBu'),0.1))
    axis equal
    hold on
    fill(obst_r*cos(0:0.01:2*pi)/d,obst_r*sin(0:0.01:2*pi)/d,'m','edgecolor','none')
    xlim([-2,10+0*max((x(:)-obst_x)/d)])
    title(sprintf('Mode %d',jj-1),'interpreter','latex','fontsize',14)
    exportgraphics(gcf,sprintf('saved_figures/compressed_modes%d.pdf',jj),'BackgroundColor','none','Resolution',400)
end






function [C,Theta] = rand_meas(p,n1,n2)
    % Creates random measurement matrix for the compression
    C = zeros(p,n1*n2);
    Theta = zeros(p,n1*n2);
    for ii=1:p
        % xmeas= zeros(n1,n2);
        % xmeas(ceil(n1*rand),ceil(n2*rand)) = 1;
        xmeas = randn(n1,n2);
        C(ii,:) = reshape(xmeas,n1*n2,1);
        Theta(ii,:) = reshape(fftshift(fft2(xmeas)),n1*n2,1)'/sqrt(n1*n2);
    end
end

function [W,LAM,K,U,V,Sinv] = standard_DMD(X,Y,r)
    % Applies exact DMD
    [U,S,V] = svd(X,'econ');
    r = min(rank(S),r);
    U = U(:,1:r); V = V(:,1:r); S = S(1:r,1:r); Sinv = diag(1./diag(S));
    K = (U')*Y*V*Sinv;
    [W,LAM] = eig(K,'vector');
end
