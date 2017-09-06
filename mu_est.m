function [F] = mu_est(miu,func)
% F = 0;
% n1=size(func,1);
F = mean(exp(miu*func'));
% for i=1:n1-1
%     F = F+  (1/n1)*exp(miu*func(i,1)');
% end
% F = log(F/n+1)+kappa*(miu*miu'); 
end




