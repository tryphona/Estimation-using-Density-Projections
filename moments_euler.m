%% Copyright to Andreas Tryphonides, 28.8.2017
%%
function [v] = moments_euler(beta,simul_x)
b = beta;
n = length(simul_x);
m1 = moment_euler_fun(simul_x,b);
m = 1/n*sum(m1)';
J= 1/n*sum(moment_euler_der(simul_x,b));
% J =1;
V = cov(m1);
v = (m'*V^-1*m);
end

