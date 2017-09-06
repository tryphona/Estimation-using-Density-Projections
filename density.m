%% Copyright to Andreas Tryphonides, 28.8.2017
%%
function [v] = density(optims,simx,simu,params,params_l_t,fmu_t,params_l_a,fmu_a,k)
b =optims(1);
n = length(simu);
gamma = 3;
mu_t = pchip(params,fmu_t,b);
lambda_t = pchip(params,params_l_t,b);
mu_a = pchip(params,fmu_a,b);
lambda_a = pchip(params,params_l_a,b);
switch k
    case 1
    a = 7;    
    v = -1/n*(sum(log((tpdf(simu-gamma-simx-0.01*simx.^2,a).*exp(mu_t*(moment_fun(simu,simx,b))+lambda_t)))));
    case 2
        
%Experiment 1        
%    a = 4;    
 %   v = -1/n*(sum(log((normpdf(simu,gamma,a)).*exp(mu_a*(moment_fun(simu,simx,b))+lambda_a))));
%Experiment 2
 %     a=4;
 %    v = -1/n*(sum(log((normpdf(simu-gamma-simx-0.01*simx.^2,gamma,a)).*exp(mu_a*(moment_fun(simu,simx,b))+lambda_a))));
%Experiment 3
       a = 7;    
       v = -1/n*(sum(log((tpdf(simu-gamma,a).*exp(mu_t*(moment_fun(simu,simx,b))+lambda_t)))));

end;

end

