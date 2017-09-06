%% Copyright to Andreas Tryphonides, 28.8.2017
%%
function [v] = density_euler(optims,simx,simu,p1,p2,beta_grid,params_l_t,fmu_t,params_l_a,fmu_a,k)
b =optims(1);
a = optims(2);

% rhoc = 0 ; rhor = 0.95; rhocr =0.05;
rhoc = 0 ; rhor = 0.95; rhocr =0.05;
rho = [rhoc rhocr; 0 rhor];
sigma_m = [0.05 a ;a 0.05];
% sigma_m = [0.05+a 0 ; 0 0.05+a];
[~,P] = chol(sigma_m);
if P>1
    v = inf;
    return
else
sigma = [0.05 0.002 ;0.002 0.05];

n = size(simu,2);

mu_t = interpn(beta_grid,fmu_t,b);
lambda_t = interpn(beta_grid,params_l_t,b);

mu_a = interpn(p1,p2,fmu_a,b,a,'cubic');
lambda_a = interpn(p1,p2,params_l_a,b,a,'cubic');
% mu_a = pchip(params,fmu_a,b);
% lambda_a = pchip(params,params_l_a,b);
% v = 1/n*(sum(log(normpdf(simul_u,simul_u,0.5).*exp(mu*(simul_u.^-b-2*b*simx.*simul_u)+lambda))));
switch k
    case 1
%     v = -1/n*(sum(log((mvtpdf((simx(:,2:end)-rho*simx(:,1:end-1))',sigma_m,7).*exp(mu_a*(moment_euler_fun(simx',b))+lambda_a)))));
%      v = -1/n*(sum(log((mvnpdf((simx(:,2:end)-rho*simx(:,1:end-1))',[0 0],sigma_m).*exp(mu_a*(moment_euler_fun(simx',b))+lambda_a)))));
 

a=0.01; sigma_m = [0.05 0 ; 0 0.05]; 
 mu_a = interpn(p1,p2,fmu_a,b,a,'cubic');
lambda_a = interpn(p1,p2,params_l_a,b,a,'cubic');
v = -1/n*(sum(log((mvnpdf((simx(:,2:end)-rho*simx(:,1:end-1))',[0 0],sigma_m).*exp(mu_a*(moment_euler_fun(simx',b))+lambda_a)))));

    case 2
%    v = -1/n*(sum(log((mvnpdf((simx(:,2:end)-rho*simx(:,1:end-1))',[0 0],sigma)).*exp(mu_t*(moment_euler_fun(simx',b))+lambda_t))));
%v = -1/n*(sum(log((mvtpdf((simx(:,2:end)-rho*simx(:,1:end-1))',sigma,7)).*exp(mu_t*(moment_euler_fun(simx',b))+lambda_t))));
sigma_m = [0.05 a ; a 0.05]; 
v = -1/n*(sum(log((mvnpdf((simx(:,2:end)-rho*simx(:,1:end-1))',[0 0],sigma_m).*exp(mu_a*(moment_euler_fun(simx',b))+lambda_a)))));
end;
end
if v==inf
 disp('inf')
end
end

