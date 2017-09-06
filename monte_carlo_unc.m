%% Copyright to Andreas Tryphonides, 28.8.2017
%%
clear all
options = optimset('Display','off','TolFun',1e-6,'MaxIter',10000,'TolX',1e-6,'MaxFunEvals',10000);
options=optimoptions('fminunc','Algorithm','quasi-newton','OptimalityTolerance',1e-6,'StepTolerance',1e-6,'MaxIterations',10000,'Display','off');

addpath('U:\Documents\research\research\first chapter\first chapter\first chapter matlab files\first chapter matlab files\paper_pics');
ns=50000; %simulation size for projection
n=200; %sample size
mc =1000; % monte carlo replications
% parameters 
gamma = 3; delta=1;
ss =10;
kk=1;
n0=20;
%  shocks_u = rand(ns*1000,1)-0.5;
a1=2; a2=5;

A=RandStream.create('mrg32k3a','NumStreams',mc);

%% simulating different variables for different MC cases: comment out the corresponding block here and also in density.m
% 1 
% shocks_u = tinv(rand(A,ns,mc),7);
% simul_x = gaminv(rand(A,ns,mc),2,5)-10;
% simul_u = gamma+shocks_u+delta*simul_x+0.01*simul_x.^2;
% shocks_n = norminv(rand(A,ns,mc),0,4);
% simul_n =gamma+shocks_n;

% % 2 
% shocks_u = tinv(rand(A,ns,mc),7);
% simul_x = gaminv(rand(A,ns,mc),2,5)-10;
% simul_u = gamma+shocks_u+delta*simul_x+0.01*simul_x.^2;
% shocks_n = norminv(rand(A,ns,mc),0,4);
% simul_n =gamma+shocks_n+delta*simul_x+0.01*simul_x.^2;

% 3 
shocks_u = tinv(rand(A,ns,mc),7);
simul_x = gaminv(rand(A,ns,mc),2,5)-10;
simul_u = gamma+shocks_u+delta*simul_x+0.01*simul_x.^2;
shocks_n = gaminv(rand(A,ns,mc),2,5)-10;
simul_n = gamma+shocks_n;
%%
sample = logspace(log10(n0),log10(n+n0),ss);
beta = fminsearch(@(beta) moments(beta,simul_x(:,1),simul_u(:,1)),-0.01);
%beta = 0.6978;
beta = 0.0173; % t distr
%  beta = -0.3366; % uniform
 
MSE = zeros(kk+1,ss,mc);

betas = -0.1:0.01:0.1;
parfor ppp = 1:length(betas)
    funct(ppp)= moments(betas(ppp),simul_x(:,1),simul_u(:,1));
end

% simul_n = simul_u;

%% monte carlo for GMM
for j = 1:length(sample)
     s=sample(j);
parfor i=1:mc
    
 beta_hats(i) = fminsearch(@(beta) moments(beta,simul_x(1:round(s,-1),i),simul_u(1:round(s,-1),i)),0.5);
  MSE(3,j,i)=  (beta_hats(i) - beta)^2;
 
end
end

params = [linspace(min(beta_hats)*0.6,max(beta_hats)*1.5,50)'];

%% Computing mu and lambda for different parameter values
sp = size(params);
fmu_t = NaN(sp);
fmu_a=NaN(sp);
params_l_t = NaN(sp);
params_l_a = NaN(sp);


tic
parfor kkk=1:length(params)
    func = moment_fun(simul_u(:,1),simul_x(:,1),params(kkk));
    [fmu_t(kkk,:), params_l_t(kkk,:)] = fminsearch(@(miu) mu_est(miu,func),-0.1);
    params_l_t(kkk,:)=1-log(params_l_t(kkk,:));
    
    funcs = moment_fun(simul_n(:,1),simul_x(:,1),params(kkk));
    [fmu_a(kkk,:),params_l_a(kkk,:)] = fminsearch(@(miu) mu_est(miu,funcs),-2);
    params_l_a(kkk,:)=1-log(params_l_a(kkk,:));

end
 toc     


% ls = length(params);
% i=1;
% while i<=ls
%     if abs(fmu_t(i))>50
%         params = params([1:i-1,i+1:length(params)]);
%         fmu_t = fmu_t([1:i-1,(i+1):length(fmu_t)]);
%         fmu_a = fmu_a([1:i-1,(i+1):length(fmu_a)]);
%         params_l_t = params_l_t([1:i-1,i+1:length(params_l_t)]);
%         params_l_a = params_l_a([1:i-1,i+1:length(params_l_a)]);
%         ls=ls-1;
%         i=i-1;
%     end
%     i=i+1;
% end



beta_aux = zeros(kk+1,ss,mc);


tic
% simul_n = zeros(kk+1,ns,mc);

%% Monte Carlo 
for k=1:kk+1
%      simul_n(k,:,:) = gamma+shocks_n+delta*simul_x+0.01*simul_x.^2;
%    simul_n(k,:,:) = gamma+shocks_n;    
%    options = optimset('Display','iter','TolFun',1e-6,'MaxIter',10000,'TolX',1e-6,'MaxFunEvals',10000);
  
for i = 1:length(sample)
 s=sample(i)
 simp = simul_u(1:round(s,-1),:);
 simxp= simul_x(1:round(s,-1),:);
parfor mcc = 1:mc 
     sim = simp(:,mcc);
     simx= simxp(:,mcc);
%      [beta_hat_2, val2] = fminsearch(@(beta) moments_s(beta,simx,sim,squeeze(simul_n(k,:,mcc))',simul_x(:,mcc),k,params,fmu),0.5,options);
   [beta_hat_2 val] = fminunc(@(optim) density(optim,simx,sim,params,params_l_t,fmu_t,params_l_a,fmu_a,k),[mean(params)],options);
    MSE(k,i,mcc)= (beta_hat_2(1)-beta).^2;
%    beta_aux(k,i,mcc)= beta_hat_2(2);

 end
end
end
toc
%% plotting 
text('Interpreter','latex')
h1 = figure;
Col = {'black','g','r','m'};
subplot(1,2,1)
for k=1:kk+2
    plot([round(sample(1:end),-1)],mean(squeeze(MSE(k,1:end,:)),2)','Color',eval(['Col{',num2str(k),'}']),...
        'MarkerSize',10,'LineWidth',2, 'LineStyle','--');
   
    hold all
end;
% plot([floor(sample)],mean(squeeze(correction(2,:,:)).^2,2),'Color','r','MarkerSize',10,'LineWidth',2, 'LineStyle','--');
title 'MSE of using mispecified density'
legend ('True','Base 1','CU-GMM');
xlabel 'Sample size'
ylabel 'MSE'

clear legend
subplot(1,2,2)
[d, x]=ksdensity(mean(simul_u,2),'kernel','epanechnikov');
plot(x,d,'black',...
        'MarkerSize',10,'LineWidth',2, 'LineStyle','--');

% simul_n =gamma+mean(beta_aux(1,end,:))*simul_x+norminv(rand(A,ns,mc),0,mean(beta_aux(1,end,:)));



hold all
for k=2:kk+1
[d, x]=ksdensity(mean(squeeze(simul_n),2));
plot(x,d,eval(['Col{',num2str(k),'}']),...
        'MarkerSize',10,'LineWidth',2, 'LineStyle','--');
% legend('Base')
hold all
end
legend ('True','Base');
title 'Densities for u_{i}'

save MSE_3 MSE
 
 saveas(h1,'Pictures/MSE_3.fig')
 saveas(h1,'Pictures/MSE_3.png')
