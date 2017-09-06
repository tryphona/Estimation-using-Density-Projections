%% Copyright to Andreas Tryphonides, 28.8.2017
%%
function [v] = moments(beta,simul_x,simul_u)
b = beta;
n = length(simul_u);
m1 = moment_fun(simul_u,simul_x,b);
m = 1/n*sum(m1)';
% J= 1/n*sum(moment_der(simul_u,simul_x,b));
V = cov(m1);
v = (m'*V^-1*m);
v = (m'*1^-1*m);
end

