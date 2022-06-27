function ret=Mutation(pml,pmu,fitness,lenchrom,chrom,sizepop,pop,bound)
%% Mutation operation
% pmutation             input  : mutation probability
% lenchrom              input  : the number of parameters
% chrom                 input  
% sizepop               input  : population size
% pop                   input  
% ret                   output : parameter set after mutation

fitmin=min(fitness);    % Find out the best fitness(for adaptive mutation probability)
fitmax=max(fitness);    % Find out the worst fitness
for i=1:sizepop  
    % The mutation operation is determined by mutation probability
    pick=rand;
    % Calculate adaptive mutation probability
    if fitness(i)>mean(fitness)         
        pmutation=pml;
    else
        pmutation=pml+(pmu-pml)*(fitness(i)-fitmin)/(fitmax-fitmin);
    end
    if pick>pmutation
        continue;
    end
    flag=0;
    while flag==0
        % Select the mutation position randomly
        pick=rand;
        while pick==0
            pick=rand;
        end
        pos=ceil(pick*sum(lenchrom));  
        v=chrom(i,pos);
        v1=v-bound(pos,1);
        v2=bound(pos,2)-v;
        pick=rand; 
        if pos>3
            if pick>0.5
                delta=v2*(1-pick^((1-pop(1)/pop(2))^2));
                chrom(i,pos)=v+delta;
            else
                delta=v1*(1-pick^((1-pop(1)/pop(2))^2));
                chrom(i,pos)=v-delta;
            end  
        else
            if pick>0.5
                delta=v2*(1-pick^((1-pop(1)/pop(2))^2));
                chrom(i,pos)=ceil(v+delta);
            else
                delta=v1*(1-pick^((1-pop(1)/pop(2))^2));
                chrom(i,pos)=ceil(v-delta);
            end   
        end
        % Check the availability
        flag=test(bound,chrom(i,:));     
    end
end
ret=chrom;
