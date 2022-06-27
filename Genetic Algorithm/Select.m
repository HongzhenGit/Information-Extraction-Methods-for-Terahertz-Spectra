function ret=Select(individuals,fitness,sizepop)
%% Select operation to get the next generation population
% individuals input  
% fitness     input  
% sizepop     input  : Population Size
% ret         output 

fitness= 1./(fitness);
sumfitness=sum(fitness);
sumf=fitness./sumfitness;
index=[];
for i=1:sizepop  
    % Perform selection operation for sizepop times
    pick=rand;
    while pick==0
        pick=rand;
    end
    % During the selection, one individual might be selected for multiple
    % times, which will ruduce the variety of population
    for j=1:sizepop
        pick=pick-sumf(j);
        if pick<0
            index=[index j];
            break;  
        end         
    end
end
individuals=individuals(index,:);
fitness=fitness(index);
ret=individuals;