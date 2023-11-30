clear
close all

addpath(genpath('./data_sets'))
addpath(genpath('./routines'))
addpath(genpath('./saved_data_from_runs'))

%% UNCOMMENT TO RUN CODE
% % load('Cylinder_wake_data.mat')
% % DATA = DATA(1:160000,:);
% % DATA = (DATA- mean(DATA,2) )./std(DATA,[],2);%
% % %% Set the parameters
% % Mvec = 24*(5:41);
% % sigma = 0.4;
% % Ns = 1;%100;
% % r = 11; % rank
% % 
% % %% Compute Error
% % rng(1)
% % ct = 1;
% % X = DATA(:,1:max(Mvec));
% % Y = DATA(:,2:(max(Mvec)+1));
% % 
% % 
% % [~,LAM] =tlsDMD(X,Y,15);
% % E = imag(log(LAM))*1i;
% % [~,I] = sort(abs(E),'ascend');
% % E = E(I(1:r));
% % Eval = cell(length(Mvec),Ns,5);
% % DIST = zeros(length(Mvec),Ns,6);
% % 
% % 
% % for M = Mvec
% %     X = DATA(:,1:M);
% %     Y = DATA(:,2:(M+1)); 
% % 
% %     for jj = 1:Ns
% %         jj
% %         ct
% %         Nr = sigma*randn(size(DATA,1),M+1);
% %         X = DATA(:,1:M) + Nr(:,1:end-1);
% %         Y = DATA(:,2:(M+1)) + Nr(:,2:end);
% % 
% %         [~,LAM] = exactDMD(X,Y,r);
% %         LAM = log(LAM);
% %         [~,I] = sort(abs(LAM),'ascend');
% %         LAM = LAM(I(1:r));
% %         Eval{ct,jj,1} = LAM;
% %         DIST(ct,jj,1) = VecDist(E,LAM);
% % 
% %         [~,LAM] = fbDMD(X,Y,r);
% %         LAM = log(LAM);
% %         Eval{ct,jj,2} = LAM;
% %         DIST(ct,jj,2) = VecDist(E,LAM);
% % 
% %         [~,LAM] = tlsDMD(X,Y,r);        
% %         LAM = log(LAM);
% %         Eval{ct,jj,3} = LAM;
% %         DIST(ct,jj,3) = VecDist(E,LAM);
% % 
% %         [~,LAM] = optdmd([X,Y(:,end)],(0:M)/10,r,2,varpro_opts('ifprint',0),1i*imag(LAM)*10);
% %         Eval{ct,jj,4} = LAM/10;
% %         DIST(ct,jj,4) = VecDist(E,LAM/10);
% %     end
% %     ct = ct+1
% % end



%%
load('Noisy_Cylinder_Results.mat')
figure

mm = mean(DIST(:,:,1),2);
st = std(DIST(:,:,1),[],2);
h=loglog(Mvec(:),mm,'k');
hold on
errorbar(Mvec(:),mm,st,'.-k','linewidth',1,'markersize',12)
h. HandleVisibility='off';

mm = mean(DIST(:,:,2),2);
st = std(DIST(:,:,2),[],2);
errorbar(Mvec(:),mm,st,'.-g','linewidth',1,'markersize',12)

mm = mean(DIST(:,:,3),2);
st = std(DIST(:,:,3),[],2);
errorbar(Mvec(:),mm,st,'.-r','linewidth',1,'markersize',12)

mm = mean(DIST(:,:,4),2);
st = std(DIST(:,:,4),[],2);
errorbar(Mvec(:),mm,st,'.-b','linewidth',1,'markersize',12)
xlim([100,1000])
ylim([10^(-5),2])

title('Mean Relative Eigenvalue Error','interpreter','latex','fontsize',18)
xlabel('$M$','interpreter','latex','fontsize',18)
grid on
ax=gca; ax.FontSize=18;

legend({'exactDMD','fbDMD','tlsDMD','optDMD'},'interpreter','latex','fontsize',12,'location','southwest')
exportgraphics(gcf,'cylinder_denoised.pdf','ContentType','vector','BackgroundColor','none')

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

% a = sqrt(norm(P-Q)^2/length(P));
a = norm(P(:)-Q(:))/norm(P(:));


end
