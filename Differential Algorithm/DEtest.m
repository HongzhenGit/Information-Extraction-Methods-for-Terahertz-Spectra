function flag=DEtest(bound,ret)
%% Availability check
% bound: input, the bound of parameters
% ret:   input, a parameter vector
flag=1;
n=length(ret);
for i=1:n
    if ret(i)<bound(i,1)||ret(i)>bound(i,2)
        flag=0;
    end
end
