function pop=DECode(NP,D,bound)
%% Generate parameter vectors within a specific bound
% NP:    input, the population size
% D:     input, the dimension of parameter vector
% bound: input, the bound of parameters
poptem=zeros(NP,D);
for i=1:NP
    flag=0;
    while flag==0
        for j=1:D
            poptem(i,j)=bound(j,1)+rand*(bound(j,2)-bound(j,1));
        end
        flag=DEtest(bound,poptem(i,:));
    end    
end
pop=poptem;