%% Copyright to Andreas Tryphonides, 28.8.2017
%%
function [der] = moment_euler_der(simul_x,params)
 c_1 = simul_x(1:end-1,1);
der =-c_1/params^2;
end

