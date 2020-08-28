%% EXTACT CELL DATA for FITTING
cellType = 'LSC';% 'LSC'; % 'SCC'
C = load(['pathFileNew_',cellType,'\result_2005_2017.mat']);
RESULT = C.RESULT;
clear C;
% FILTER OUT PROBLEMATIC CELLPATH
func_trim = @(A)find(A>prctile(A,0.25) & A<prctile(A,99.75));
    
RESULT = RESULT(func_trim(RESULT.vel),:);
RESULT = RESULT(RESULT.durT > 1 | ~RESULT.newstart,:);
ii = find(RESULT.angle>180);
if ~isempty(ii)
    RESULT.angle(ii) = RESULT.angle(ii)-360;
end
cellPath = RESULT;

%% compute Kendall tau for all statistics
A = [cellPath.Axismin,cellPath.Axismaj,cellPath.ecc,...
    cellPath.meanR,cellPath.maxR,cellPath.mmaxR,...
    cellPath.vel,cellPath.angle,...
    cellPath.peakX,cellPath.peakY,...
    cellPath.cellT-cellPath.stormT,...
    cellPath.durT];%,totalTimes_storm];

% kendallTau = corr(A,'type','Kendall');
[r,p] = corr(smaj,smin,'type','Kendall');
figure;
scatterhist(smaj,smin);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------3D Copula--------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute 3D copula for mmaxR, smaj, pathLength

A = [cellPath.mmaxR,...
    cellPath.Axismaj,...
    cellPath.durT];% 1+1./(1-cellPath.ecc),...

[RHO,PVAL] = corr(A,'type','Kendall');

%% ------------------------------------------------------------------------%
% Have discarded the SUSPECT path.
% copula fit for mmaxR - Axismaj - pathLength

% mmaxR:       gamma
% Axismaj:     logNormal
% pathLength:  logNormal

% copula type: t-copula
% -------------------------------------------------------------------------%
% For 4D copula: load copula_4D_param.mat param; % for mmaxR, Axismaj, 1/1-ecc, pathLength

param = [1.81776244814184,3.99936970775806;
    2.19381480694163,0.456145826153028;
    0.813023327672115,0.829317405252584];

X = cellPath.mmaxR-35;
shape_X = param(1,1);
scale_X = param(1,2);

Y = cellPath.Axismaj;
mu_Y = param(2,1);
sigma_Y = param(2,2);

Z = cellPath.durT;
aa = randl(logncdf(Z+1,mu_Y,sigma_Y)-logncdf(Z,mu_Y,sigma_Y),size(cellPath.durT));
Z = Z+aa;% add random term for cell path length;
mu_Z = param(3,1);
sigma_Z = param(3,2);


%% TRASFORM: Distribution 2 Uniform
% option1: icdf
% u = logncdf(X,mu_X,sigma_X);
u = gamcdf(X,shape_X,scale_X);
v = logncdf(Y,mu_Y,sigma_Y);
w = ksdensity(Z,Z,'function','cdf');

% %% option2: kernell
% u = ksdensity(X,X,'function','cdf');
% v = ksdensity(Y,Y,'function','cdf');
% w = ksdensity(Z,Z,'function','cdf');

save([dataPath,'icdf_3D_copula_0304.mat'],'X','Y','Z','u','v','w','param','A');

% %% CHECKING: check whether normal
% figure;
% scatterhist(v,w)
% axis([0 1 0 1]);
% xlabel('U(***)');
% ylabel('U(pathLength)');
% title('U(***)-U(pathLength)');
% set(get(gca,'children'),'marker','.','markersize',1);

%% COPULA FITTING
rng default  % For reproducibility
[Rho_fit,Nu_fit] = copulafit('t',[u v w]);

%% COPULA GENERATION
r = copularnd('t',Rho_fit,Nu_fit,length(u));
% r = copularnd('Clayton',paramhat,length(u));
% r = copularnd('Frank',paramhat,length(u));

u1 = r(:,1);
v1 = r(:,2);
w1 = r(:,3);

% Tranformation
% X Y Z
% X1 Y1 Z1
X1 = gaminv(u1,shape_X,scale_X);
Y1 = exp(norminv(v1) .* sigma_Y + mu_Y);
Z1 = exp(norminv(w1) .* sigma_Z + mu_Z);


%% VALIDATION: PLOT
%% mmaxR vs pathLength
ii = find(X1 < 35 & Z1 < 50);

ZZ1 = floor(Z1(ii));
ZZ1(ZZ1 == 0) = 1;
figure;
scatterhist(Y1(ii),ZZ1);
set(get(gca,'children'),'marker','.');
xlabel('sim(mmaxR)');
ylabel('sim(pathLength)');
title('sim(mmaxR)-sim(pathLength)');
set(get(gca,'children'),'marker','.');
set(gca,'xlim',[0 60],'ylim',[0 50]);

Z_ori = cellPath.durT;
ii = find(X < 35 & Z < 50);
figure;
scatterhist(Y(ii),floor(Z_ori(ii)));
set(get(gca,'children'),'marker','.');
xlabel('obs(mmaxR)');
ylabel('obs(pathLength)');
title('obs(mmaxR)-obs(pathLength)');
set(get(gca,'children'),'marker','.');
set(gca,'xlim',[0 60],'ylim',[0 50]);

%% mmaxR vs AxisMaj
ii = find(X1 < 65);% & Y1 < 50);

% ii = ii(1:length(X));
figure;
scatterhist(X1(ii),Y1(ii));
set(get(gca,'children'),'marker','.');
xlabel('sim(mmaxR)');
ylabel('sim(Axismaj)');
title('sim(mmaxR)-sim(Axismaj)');
set(get(gca,'children'),'marker','.');
set(gca,'xlim',[0 40],'ylim',[0 60]);

ii = find(X < 35);% & F < 50);
figure;
scatterhist(X(ii),Y(ii));
set(get(gca,'children'),'marker','.');
xlabel('obs(mmaxR)');
ylabel('obs(Axismaj)');
title('obs(mmaxR)-obs(Axismaj)');
set(get(gca,'children'),'marker','.');
set(gca,'xlim',[0 40],'ylim',[0 60]);

%%
ii = find(Z1 < 50);
ZZ1 = floor(Z1(ii));
ZZ1(ZZ1 == 0) = 1;
A = [X1(ii),Y1(ii),Z1(ii),ZZ1];

%% VALIDATION: STATISTICS

[Kendall_fit,PVAL_fit] = corr(A,'type','Kendall');
[Pearson_fit,PVAL_fit2] = corr(A,'type','Pearson');
