%% Copyright to Andreas Tryphonides, 28.8.2017
%%
function [der] = moment_der(simul_u,simul_x,params)
   der = (-params)*(squeeze(simul_u.^(-params-1))-2*squeeze(simul_u).*simul_x);
end

