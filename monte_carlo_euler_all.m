%% Copyright to Andreas Tryphonides, 28.8.2017
%%
% clf
options = optimset('Display','off','TolFun',1e-4,'MaxIter',10000,'TolX',1e-4,'MaxFunEvals',10000);
ns=500;
n=200;
mc =500;
ss =10;
kk=1;
n0=10;
sample = logspace(log10(n0),log10(n+n0),ss);

%% Simulating variables for different MC cases (different restrictions on the base density) 
rhoc = 0 ; rhor = 0.95; rhocr = 0.05;
rho = [rhoc rhocr; 0 rhor];
sigma = [0.05 0.002 ;0.002 0.05];

A=RandStream.create('mrg32k3a','NumStreams',mc);

rng(2000);
ls = 10;
sig12_grid = linspace(0.0002,0.01,ls);
delta = 0.2;

rho12_grid = linspace(0.005,0.15,ls);
rho21_grid = linspace(-0.001,0.001,ls);

rho12 = rhocr*(1-delta./sqrt(linspace(2,ns,ls)));
rho21 = delta./sqrt(linspace(2,ns,ls));
rho22 = 1-rho12;
rho11 = rhoc;

simul_x = zeros(2,ns,mc);
simul_u = zeros(2,ns,ls,ls,ls,mc);
% 
parfor j=1:mc
    simul_xs = ones(2,ns)*0.001;
    simul_us = ones(2,ns,ls,ls,ls)*0.001;
    j
for i=2:ns
%     c = squeeze(simul_xs(1,i-1));
%     r = squeeze(simul_xs(2,i-1));
    for k = 1:ls
        for l=1:ls
            for m = 1:ls            
       sigma_m = [0.05 sig12_grid(k) ;sig12_grid(k) 0.05]; 
 %       errs_us = mvtrnd(sigma_m,7)';
      errs_us = mvnrnd([0 0]',sigma_m)';
     rho_u=[rhoc 0;0 rhor];
     rho_u=[rho11 rho12_grid(l) ; rho21_grid(m) 1-rho12_grid(l)];         
     simul_us(:,i,k,l,m) = rho_u*simul_us(:,i-1,k,l,m)+errs_us;
     
        errs_xs = mvnrnd([0 0]',sigma)';
%     errs_xs = mvtrnd(sigma,15)';
    simul_xs(:,i) = rho*simul_xs(:,i-1)+errs_xs;
     
     
     
            end
        end
    end

    
end
simul_u(:,:,:,:,:,j) = simul_us; 
simul_x(:,:,j) = simul_xs;
end
  

save('simul_mc_all.mat','simul_u','simul_x','-v7.3')
%%
load simul_mc_all



% load simul_mc_rho
%load simul_mc

betat = 0.7447;% true value when normal distribution is used
%betat = 9.3750e-04;% true value when t distribution is used

% simul_x_for = simul_x(:,2:end,:,:,:,:);
% simul_x_lag = simul_x(:,1:end-1,:,:,:,:);

%% Monte Carlo for GMM
MSE = zeros(kk+1,ss,mc);
for ss=1:length(sample)
    s=sample(ss);
parfor i=1:mc
%   c = simul_x_for(1,2:round(s,-1),i);
%   c_1 = simul_x_lag(1,1:round(s,-1)-1,i);
%   r = simul_x_for(2,2:round(s,-1),i);
 beta_hats(i) = fminsearch(@(beta) moments_euler(beta,squeeze(simul_x(:,1:round(s,-1),i))'),0.6);
 MSE(3,ss,i)= (beta_hats(i)-betat)^2;

end
end
%%

beta_grid = [linspace(min(min(beta_hats))-0.1,1,50)'];
% params = [beta_grid sig12_grid'];
%[params1,params2,params3,params4] = meshgrid(beta_grid,sig12_grid,rho12_grid,rho21_grid);

sp = [length(beta_grid),ls,ls,ls];
fmu_t = NaN([sp(1),mc]);
fmu_a= zeros([sp mc]);
%params_l_t = NaN([sp(1),mc]);
 params_l_a = zeros([sp mc]);
% 
% % 
% tic

%% Conputing mu and lambda  
for kkk=1:length(beta_grid)
    parfor mcc=1:mc
%     func = moment_euler_fun(simul_x(:,:,mcc)',params1(kkk,1));
%     [a, b] = fminsearch(@(miu) mu_est(miu,func),0.01);
%     fmu_t(kkk,mcc)=a;
%     params_l_t(kkk,mcc)=1-log(b);    
          
    for iii = 1:ls
        for lll =1:ls
            for mmm = 1:ls
    funcs = moment_euler_fun(squeeze(simul_u(:,:,iii,lll,mmm,mcc))',beta_grid(kkk));
    [a,b] = fminsearch(@(miu) mu_est(miu,funcs),-0.01,options);
    fmu_a(kkk,iii,lll,mmm,mcc)=a;
    params_l_a(kkk,iii,lll,mmm,mcc)=1-log(b);
            end
        end
    end
    end
end
 toc 
 %% Monte Carlo

beta_aux = zeros(kk+1,ss,mc);


tic
% simul_n = zeros(kk+1,ns,mc);
[p1,p2,p3,p4] = ndgrid(beta_grid,sig12_grid,rho12_grid,rho21_grid );
b1 = ndgrid(beta_grid);
kk = 1;
for k=1:kk+1
for i = 1:length(sample)
 s=sample(i)
 simp = squeeze(simul_u(:,1:round(s,-1),:,:,:,:));
 simxp= squeeze(simul_x(:,1:round(s,-1),:));
tic
 parfor mcc = 1:mc  
    
    sim = squeeze(simp(:,:,:,:,:,mcc));
     simx= squeeze(simxp(:,:,mcc));
  
   [beta_hat_2, val] = fminsearch(@(optim) density_euler_all(optim,simx,sim,p1,p2,p3,p4,b1,[],[],params_l_a(:,:,:,:,mcc),fmu_a(:,:,:,:,mcc),k),[0.7 0.001 0.01 0.001],options);
    MSE(k,i,mcc)= (beta_hat_2(1)-betat)^2;
 end
end
end

h1 = figure;
Col = {'black','g','r','m'};
% subplot(1,2,1)
for k=1:kk+2
    plot(round(sample(1:end),-1),mean(squeeze(MSE(k,1:end,:)),2)','Color',eval(['Col{',num2str(k),'}']),...
        'MarkerSize',10,'LineWidth',2, 'LineStyle','-.');
   
    hold all
end;
title 'MSE of using mispecified density'
% legend ('Base-True Estimated (\sigma_{CR})','Base-True Model','CU-GMM');
legend ('Base-Restrictions on (\sigma_{CR})','Base-Estimated True Model','CU-GMM');
xlabel 'Sample size'
ylabel 'MSE'


save MSE_euler_all MSE
saveas(h1,'MSE_euler_all.fig')

