%% GA For Parameter Estimation
clc 
clear all; 
close all; 
warning off
%% Initializing
popsize=100;              % Population Size
lenchrom=5;               % The number of parameters

pcl=0.1;                  % Set the lower bound of Cross Probability
pcu=0.6;                  % Set the upper bound of Cross Probability
pml=0.1;
pmu=0.3;
% pc=0.7;                 % Fixed cross probablity
% pm=0.5;                 % Fixed mutation probablity

maxgen=500;               % The maximum number of iterations

% Population features
popmax=3000;
popmin=0;
popmax1=2000;
popmin1=0;
popmax2=1;
popmin2=-1;

% The scale of parameters(A proper range of parameters could fast the convergence speed)
bound=[popmin popmax;popmin1 popmax1;popmin1 popmax1;popmin2 popmax2;popmin2 popmax2];  

% Read the refrence/measurement signals of samples
MeaData=textread('GA_measured.txt');     
ts=MeaData(:,1);
measure=MeaData(:,2);
% measure=eemd(measure,0.2,150);
N=length(measure);
RefData=textread('GA_reference.txt');                   
reference=RefData(:,2);
reference=sigshift(reference,9000,N);

% Plot the signals in frequency domain(after FFT)
% ddt=ts(2)-ts(1);
% n=0:N-1;
% f=n/(N*ddt);
% fmeasure=fft(measure);
% freference=fft(reference);
% figure
% plot(f(1:200),abs(fmeasure(1:200)),'b-','LineWidth',2);
% hold on
% plot(f(1:200),abs(freference(1:200)),'r-','LineWidth',2);

% Re-construct the signal with EMD method
% reference=eemd(reference,0.2,150);
% sum=reference(:,5);
% for i=6:11
%     sum=sum+reference(:,i);
% end
% sum1=measure(:,5);
% for i=6:11
%     sum1=sum1+measure(:,i);
% end
% reference=sum;
% measure=sum1;

%% Initiate the population
for i=1:popsize
    % Generate a population randomly
    GApop(i,:)=Code(lenchrom,bound);       
    % Calculate the fitness of each individual
    fitness(i)=fun(GApop(i,:),reference,measure,N);            
end

% Find out he best parameter vector
[bestfitness,bestindex]=min(fitness);
zbest=GApop(bestindex,:);   % Population Best
gbest=GApop;                % Individual Best
fitnessgbest=fitness;       % Fitness of the Population Best
fitnesszbest=bestfitness;   % Fitness of the Individual Best

%% Run the iterations
for i=1:maxgen
        i
        
        % Population update
        GApop=Select(GApop,fitness,popsize);
        % Cross
        GApop=Cross(pcl,pcu,fitness,lenchrom,GApop,popsize,bound);
        % Mutation
        GApop=Mutation(pml,pmu,fitness,lenchrom,GApop,popsize,[i maxgen],bound);

        pop=GApop;
        
      for j=1:popsize
        % Fitness caculation
        fitness(j)=fun(pop(j,:),reference,measure,N);
        % Individual Best Update
        if fitness(j) < fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
        
        % Population Best Update
        if fitness(j) < fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j);
        end        
     end
    
    yy(i)=fitnesszbest;     
end

%% Print the results
disp '**************best parameter set*****************'
zbest
% Calculate the desired information
% Flight time, reflection coeff, reflection index and thickness d.
dt=abs(zbest(3)-zbest(2))*(ts(2)-ts(1));
r=max(zbest(4),zbest(5));
n=(1+r)/(1-r);
d=((3*10^8*dt*10^(-12))/(2*n))*10^6;

results=[dt,r,n,d];

disp '*****************   results   *******************'
results
% Final fitted results
cankaoshift=sigshift(reference,zbest(1),N);
sig1=zbest(4)*sigshift(cankaoshift,zbest(2),N);
sig2=zbest(5)*sigshift(cankaoshift,zbest(3),N);
% sig3=zzbest(6)*sigshift(cankaoshift,max(zzbest(2),zzbest(3))+abs(zzbest(2)-zzbest(3)),N);
sig=sig1+sig2;
%%
plot(yy,'linewidth',2);
title(['Fiteness Curve' 'Iteration Stopped£½' num2str(maxgen)]);
xlabel('Number of iterations');ylabel('Fitness');
grid on
figure
plot(ts,sig1,'r-','LineWidth',2)
hold on
plot(ts,sig2,'b-','LineWidth',2)
% hold on
% plot(ts,sig3,'g-','LineWidth',2)
hold on
plot(ts,sig,'k-','LineWidth',2)
xlabel('Time/ps');
ylabel('Amplitude');
title('Fitted signal and its components');
legend('1st Refelection','2nd Reflection','Fitted Result');
figure
plot(ts,sig,'r-','LineWidth',2)
hold on
plot(ts,measure,'b-','LineWidth',2)
xlabel('Time/ps');
ylabel('Amplitude');
title('Fitted Signal VS. Actual Signal');
legend('Fitted Signal','Actual Signal');


