%% Copyright to Andreas Tryphonides, 28.8.2017
%%
function [v] = density_euler_all(optims,simx,simu,p1,p2,p3,p4,beta_grid,params_l_t,fmu_t,params_l_a,fmu_a,k)
b =optims(1);
a1 = optims(2);
a2 = optims(3);
a3 = optims(4);


% rhoc = 0 ; rhor = 0.95; rhocr =0.05;
rhoc = 0 ; rhocr =a2;rhorc=a3;
rho = [rhoc rhocr; rhorc 1-rhocr];
sigma = [0.05 0.002 ;0.002 0.05];
sigma_m = [0.05 a1 ;a1 0.05];
[~,P] = chol(sigma_m);
if P>1
    v = inf;
    return
else

n = size(simu,2);
% mu_t = pchip(params,fmu_t,b);

% mu_t = interpn(beta_grid,fmu_t,b);
% lambda_t = interpn(beta_grid,params_l_t,b);


% mu_a = pchip(params,fmu_a,b);
% lambda_a = pchip(params,params_l_a,b);
% v = 1/n*(sum(log(normpdf(simul_u,simul_u,0.5).*exp(mu*(simul_u.^-b-2*b*simx.*simul_u)+lambda))));
switch k
    case 1
%     v = -1/n*(sum(log((mvtpdf((simx(:,2:end)-rho*simx(:,1:end-1))',sigma_m,7).*exp(mu_a*(moment_euler_fun(simx',b))+lambda_a)))));
%      v = -1/n*(sum(log((mvnpdf((simx(:,2:end)-rho*simx(:,1:end-1))',[0 0],sigma_m).*exp(mu_a*(moment_euler_fun(simx',b))+lambda_a)))));
 
    mu_a = interpn(p1,p2,p3,p4,fmu_a,b,a1,a2,a3,'cubic');
if isnan(mu_a)==1
    v=inf;
    return
end
lambda_a = interpn(p1,p2,p3,p4,params_l_a,b,a1,a2,a3,'cubic');



 
v = -1/n*(sum(log((mvnpdf((simx(:,2:end)-rho*simx(:,1:end-1))',[0 0],sigma_m).*exp(mu_a*(moment_euler_fun(simx',b))+lambda_a)))));

case 2
   a2= 0.05+0.01/n;
   a3=0.01/n;
   
    rho = [rhoc a2; a3 1-a2];
    mu_a = interpn(p1,p2,p3,p4,fmu_a,b,a1,a2,a3,'cubic');
if isnan(mu_a)==1
    v=inf;
    return
end
lambda_a = interpn(p1,p2,p3,p4,params_l_a,b,a1,a2,a3,'cubic');

   
    
    

v = -1/n*(sum(log((mvnpdf((simx(:,2:end)-rho*simx(:,1:end-1))',[0 0],sigma_m).*exp(mu_a*(moment_euler_fun(simx',b))+lambda_a)))));
end;
end
if v==inf
 disp('inf')
end
end

