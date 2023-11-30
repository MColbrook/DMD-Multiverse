clear
close all

addpath(genpath('./data_sets'))
addpath(genpath('./routines'))

load('Cylinder_wake_data.mat') % load the data
%% Sort the physical coordinates
x=(x-obst_x)/d; y=(y-obst_y+2.5)/d; x=x-4;
[~,I]=sort(abs(x+1i*y),'ascend');

%% Collect data at a single point
N = 100;
M = 24*5;
D = DATA(I(1),1:(M+N));

%% Form the Hankel matrices
X = zeros(N,M); X(1,:) = D(1:M);
Y = zeros(N,M); Y(1,:) = D((1:M)+1);

for jj = 2:N
    X(jj,:) = D(jj:(M+jj-1));
    Y(jj,:) = D((jj+1):(M+jj));
end

%% Perform exact DMD
[~,LAM] = exactDMD(X,Y,39); 

%% Plot the singular values
S = svd(X);
figure
semilogy(S,'b.','markersize',15)
title('Singular Values','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
exportgraphics(gcf,'saved_figures/Hankel_SVD.pdf','ContentType','vector','BackgroundColor','none')

%% Plot the eigenvalues
figure
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
hold on
plot(real(LAM),imag(LAM),'b.','markersize',15)
axis equal
axis([-1.15,1.15,-1.15,1.15])
title('Hankel-DMD Eigenvalues','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
exportgraphics(gcf,'saved_figures/Hankel_eigenvalues.pdf','ContentType','vector','BackgroundColor','none')


%% Compute the error in the eigenvalues with M
Mvec = 11:100;
Er = 0*Mvec;
N=100;

M = 24*10;
ind = (1:M);
X = DATA(1:160000,ind);
Y = DATA(1:160000,ind+1);

[~,LAM0] = exactDMD(X,Y,41);
LAM0 = imag(log(LAM0))*1i;
[~,II] = sort(abs(LAM0),'ascend');
LAM0 = LAM0(II(1:11));

ct = 1;
for M = Mvec
    D = DATA(I(1),1:(M+N));
    X = zeros(N,M); X(1,:) = D(1:M);
    Y = zeros(N,M); Y(1,:) = D((1:M)+1);
    for jj = 2:N
        X(jj,:) = D(jj:(M+jj-1));
        Y(jj,:) = D((jj+1):(M+jj));
    end
    [~,LAM] = exactDMD(X,Y,39);
    LAM = log(LAM);
    [~,II] = sort(abs(LAM),'ascend');
    LAM = LAM(II(1:11));
    Er(ct) = VecDist(LAM0,LAM);
    ct=ct+1;
end

figure
semilogy(Mvec,Er,'b','linewidth',2)
xlim([10,100])

title('Relative Eigenvalue Error','interpreter','latex','fontsize',18)
xlabel('$M$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
exportgraphics(gcf,'saved_figures/Hankel_error.pdf','ContentType','vector','BackgroundColor','none')




function a = VecDist(P,Q)
%%% CURRENTLY ASSUMES LENGTH OF P AND Q ARE EQUAL %%%
P = P(:); Q = Q(:); Q = flipud(Q);

C = zeros(length(P),length(Q));
for ii = 1:length(P)
    for jj = 1:length(Q)
        C(ii,jj)=abs(P(ii)-Q(jj))^2;
    end
end

M = matchpairs(C,1000);

P = P(M(:,1));
Q = Q(M(:,2));
[~,I] = sort(abs(P),'ascend');
P = P(I(1:11));
Q = Q(I(1:11));


a = norm(P(:)-Q(:))/norm(P(:));

end
