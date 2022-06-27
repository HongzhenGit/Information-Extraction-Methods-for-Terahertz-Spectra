function flag=test(bound,code)
%% Test whether a parameter set is falling into its bound
% lenchrom   input 
% bound      input 
% code       output
flag=1;
[n,m]=size(code);

for i=1:n
    if code(i)<bound(i,1) || code(i)>bound(i,2)
        flag=0;
    end
end