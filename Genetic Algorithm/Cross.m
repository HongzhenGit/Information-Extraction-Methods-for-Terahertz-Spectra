function ret=Cross(pcl,pcu,fitness,lenchrom,chrom,sizepop,bound)
%% Cross operation
% pcorss                input  : cross probability
% lenchrom              input  : the number of parameters
% chrom                 input  
% sizepop               input  : population size
% ret                   output : the result after cross operation

fitness= 1./(fitness);       % Turns the target into a minimum optimization problem
fitmin=min(fitness);         % Find out the best fitness(for adaptive cross probability)
fitmax=max(fitness);         % Find out the worst fitness(for adaptive cross probability)

for i=1:sizepop   
    % Random perform cross operation on two parameter vectors
    pick=rand(1,2);
    while prod(pick)==0
        pick=rand(1,2);
    end
    index=ceil(pick.*sizepop);
    % Determine whether a cross operation would happen based on probability
    pick=rand;
    while pick==0
        pick=rand;
    end
    % Calculate the adaptive probability
    if fitness(i)>mean(fitness)         
        pcross=pcl+(pcu-pcl)*(fitness(i)-fitmin)/(fitmax-fitmin);
    else
        pcross=pcl;
    end
    if pick>pcross
        continue;
    end
    flag=0;
    while flag==0
        % Select the cross position randomly
        pick=rand;
        while pick==0
            pick=rand;
        end
        % The cross operation are happened at the same position
        pos=ceil(pick.*sum(lenchrom)); 
        pick=rand; 
        if pos>3
            v1=chrom(index(1),pos);
            v2=chrom(index(2),pos);
            chrom(index(1),pos)=pick*v2+(1-pick)*v1;
            chrom(index(2),pos)=pick*v1+(1-pick)*v2; 
            % Check the availability
            flag1=test(bound,chrom(index(1),:));  
            flag2=test(bound,chrom(index(2),:));  
        else
            v1=chrom(index(1),pos);
            v2=chrom(index(2),pos);
            chrom(index(1),pos)=ceil(pick*v2+(1-pick)*v1);
            chrom(index(2),pos)=ceil(pick*v1+(1-pick)*v2);
            % Check the availability
            flag1=test(bound,chrom(index(1),:));  
            flag2=test(bound,chrom(index(2),:));  
        end
        if   flag1*flag2==0
            flag=0;
        else
            flag=1;
        end    
    end
end
ret=chrom;
