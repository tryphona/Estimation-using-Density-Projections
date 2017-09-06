%% Copyright to Andreas Tryphonides, 28.8.2017
%%
function [f] = moment_euler_fun(simul_x,params)
betta=params(1);
c = simul_x(2:end,1);
  c_1 = simul_x(1:end-1,1);
   r = simul_x(2:end,2);
f= exp(c+r)-exp(c_1-log(betta));
end