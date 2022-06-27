function ret=Code(lenchrom,bound)
%% Generate parameter sets within a specific bound
% lenchrom   input : the number of parameters
% bound      input : the bound of parameters
% ret        output: final coded paramter vector

flag=0;
while flag==0
    pick=rand(1,lenchrom);
    % Random interploting
    ret=bound(:,1)'+(bound(:,2)-bound(:,1))'.*pick;
    ret=[ceil(ret(1:3)),ret(4:5)];
    % Check whether a parameter set is available 
    flag=test(bound,ret);             
end
