%% Copyright to Andreas Tryphonides, 28.8.2017
%%
function [f] = moment_fun(simul_n,simul_x,params)
  f=exp(simul_n).^(-params)-simul_x.*simul_n*params;
end