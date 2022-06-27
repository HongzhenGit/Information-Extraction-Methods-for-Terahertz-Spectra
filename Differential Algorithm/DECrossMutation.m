function ui=DECrossMutation(strategy,rot,NP,D,popold,bestmemit,bound,...
                            f1,c,m1,N,ml1,ml2,ml3,ml4,measure1,...
                            Freference1_obj,mod,paint,factor)
%% Cross and Mutation Operation 
% strategy:  input, strategy indicator
% rot:       input, a parameter vector of population
% NP:        input, population size
% popold:    input, population from last iteration
% bestmemit: input, the best vector of current iteration
% bound:     input, the bounds of parameter vector
% f1,c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod,paint,factor:
% inputs for fitness function
CRL=0.1;
CRU=0.6;

FL=0.3;
FU=0.8;
FF=0.6;

for i=1:NP                                                                 % population filled with the best member
    bm(i,:) = bestmemit;                                                   % of the last iteration
end

ind=randperm(4);                                                           % index pointer array

a1=randperm(NP);                                                           % shuffle locations of vectors 
rt=rem(rot+ind(1),NP);                                                     % rotate indices by ind(1) positions
a2=a1(rt+1);                                                               % rotate vector locations
rt=rem(rot+ind(2),NP);
a3=a2(rt+1);                
rt=rem(rot+ind(3),NP);
a4=a3(rt+1);               
rt=rem(rot+ind(4),NP);
a5=a4(rt+1);                

pm1=popold(a1,:);                                                          % shuffled population 1
pm2=popold(a2,:);                                                          % shuffled population 2
pm3=popold(a3,:);                                                          % shuffled population 3
pm4=popold(a4,:);                                                          % shuffled population 4
pm5=popold(a5,:);                                                          % shuffled population 5
if paint==1
    fitbest=Obfun(f1,bestmemit,c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod,factor);
elseif paint==2
    fitbest=SObfun(f1,bestmemit,c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
elseif paint==3
    fitbest=TObfun(f1,bestmemit,c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
else
    fitbest=FObfun(f1,bestmemit,c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
end
for i=1:NP
    if paint==1
        fitpm1(i)=Obfun(f1,pm1(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod,factor);
    elseif paint==2
        fitpm1(i)=SObfun(f1,pm1(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
    elseif paint==3
        fitpm1(i)=TObfun(f1,pm1(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
    else
        fitpm1(i)=FObfun(f1,pm1(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
    end
end
for i=1:NP
    if paint==1
        fitpm2(i)=Obfun(f1,pm2(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod,factor);
    elseif paint==2
        fitpm2(i)=SObfun(f1,pm2(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
    elseif paint==3
        fitpm2(i)=TObfun(f1,pm2(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
    else
        fitpm2(i)=FObfun(f1,pm2(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
    end
end

flag=0;

while flag==0
    if strategy==1
        fiteval=[fitpm1',fitpm2'];
        for i=1:NP
            fitmax=max(fiteval(i,:));
            fitmin=min(fiteval(i,:));
            F(i)=FL+(FU-FL)*(fitmin-fitbest)/(fitmax-fitbest);
            if ismissing(F(i))
                F(i)=0.5;
            end
            ui(i,:)=bm(i,:)+F(i)*(pm1(i,:)-pm2(i,:));
        end
        for i=1:NP
            if paint==1
                fit(i)=Obfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod,factor); 
            elseif paint==2
                fit(i)=SObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod); 
            elseif paint==3
                fit(i)=TObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
            else
                fit(i)=FObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod); 
            end
        end
        fitm=mean(fit);
        fitmax=max(fit);
        fitmin=min(fit);
        for i=1:NP
            if fit(i)<fitm
                CR(i)=CRL+(CRU-CRL)*(fitmax-fit(i))/(fitmax-fitmin);
            else
                CR(i)=CRL;
            end
        end
        CRM=[ ];
        for i=1:D
            CRM=[CRM,CR'];
        end
        mui=rand(NP,D)<CRM;
        mpo=mui < 0.5;
        ui=popold.*mpo+ui.*mui;
        for i=1:NP
             flag=DEtest(bound,ui(i,:));
             while flag==0
                 for j=1:D
                     ui(i,j)=bestmemit(1,j)+(2*rand-1)*(bound(j,2)...
                         -bestmemit(1,j));
                 end
                 flag=DEtest(bound,ui(i,:));
             end
        end
    elseif strategy==2
        for i=1:NP
            if paint==1
                fitpm3(i)=Obfun(f1,pm3(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod,factor);
            elseif paint==2
                fitpm3(i)=SObfun(f1,pm3(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
            elseif paint==3
                fitpm3(i)=TObfun(f1,pm3(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
            else
                fitpm3(i)=FObfun(f1,pm3(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
            end
        end
        fiteval=[fitpm1',fitpm2',fitpm3'];
        for i=1:NP
            [fitmax,maxindex]=max(fiteval(i,:));
            fitmedian=median(fiteval(i,:));
            medindex=find(fiteval(i,:)==median(fiteval(i,:)));
            [fitmin,minindex]=min(fiteval(i,:));
            F(i)=FL+(FU-FL)*(fitmedian-fitmin)/(fitmax-fitmin);
            if ismissing(F(i))
                F(i)=0.5;
            end
            ret=[pm1(i,:)',pm2(i,:)',pm3(i,:)'];
            ret=ret';
            ui(i,:)=ret(minindex,:)+F(i)*(ret(medindex,:)-ret(maxindex,:));
        end
        for i=1:NP
            if paint==1
                fit(i)=Obfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod,factor); 
            elseif paint==2
                fit(i)=SObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
            elseif paint==3
                fit(i)=TObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
            else
                fit(i)=FObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod); 
            end
        end
        fitm=mean(fit);
        fitmax=max(fit);
        fitmin=min(fit);
        for i=1:NP
            if fit(i)<fitm
                CR(i)=CRL+(CRU-CRL)*(fitmax-fit(i))/(fitmax-fitmin);
            else
                CR(i)=CRL;
            end
        end
        CRM=[ ];
        for i=1:D
            CRM=[CRM,CR'];
        end
        mui=rand(NP,D)<CRM;
        mpo=mui < 0.5;
        ui=popold.*mpo+ui.*mui;
        for i=1:NP
             flag=DEtest(bound,ui(i,:));
             while flag==0
                 for j=1:D
                     ui(i,j)=pm3(i,j)+(2*rand-1)*(bound(j,2)...
                         -pm3(i,j));
                 end
                 flag=DEtest(bound,ui(i,:));
             end
        end
    elseif strategy==3
        ui=popold+FF*(bm-popold)+FF*(pm1-pm2);
        for i=1:NP
            if paint==1
                fit(i)=Obfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod,factor); 
            elseif paint==2
                fit(i)=SObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod); 
            elseif paint==3
                fit(i)=TObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod); 
            else
                fit(i)=FObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod); 
            end
        end
        fitm=mean(fit);
        fitmax=max(fit);
        fitmin=min(fit);
        for i=1:NP
            if fit(i)<fitm
                CR(i)=CRL+(CRU-CRL)*(fitmax-fit(i))/(fitmax-fitmin);
            else
                CR(i)=CRL;
            end
        end
        CRM=[ ];
        for i=1:D
            CRM=[CRM,CR'];
        end
        mui=rand(NP,D)<CRM;
        mpo=mui < 0.5;
        ui=popold.*mpo+ui.*mui;
        for i=1:NP
             flag=DEtest(bound,ui(i,:));
             while flag==0
                 for j=1:D
                     ui(i,j)=popold(i,j)+(2*rand-1)*(bound(j,2)...
                         -popold(i,j));
                 end
                 flag=DEtest(bound,ui(i,:));
             end
        end
    elseif strategy==4
        ui=bm+FF*(pm1-pm2+pm3-pm4);
        for i=1:NP
            if paint==1
                fit(i)=Obfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod,factor);
            elseif paint==2
                fit(i)=SObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
            elseif paint==3
                fit(i)=TObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod); 
            else
                fit(i)=FObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod); 
            end
        end
        fitm=mean(fit);
        fitmax=max(fit);
        fitmin=min(fit);
        for i=1:NP
            if fit(i)<fitm
                CR(i)=CRL+(CRU-CRL)*(fitmax-fit(i))/(fitmax-fitmin);
            else
                CR(i)=CRL;
            end
        end
        CRM=[ ];
        for i=1:D
            CRM=[CRM,CR'];
        end
        mui=rand(NP,D)<CRM;
        mpo=mui < 0.5;
        ui=popold.*mpo+ui.*mui;
        for i=1:NP
             flag=DEtest(bound,ui(i,:));
             while flag==0
                 for j=1:D
                     ui(i,j)=bestmemit(1,j)+(2*rand-1)*(bound(j,2)...
                         -bestmemit(1,j));
                 end
                 flag=DEtest(bound,ui(i,:));
             end
        end
    else 
        ui=pm5+FF*(pm1-pm2+pm3-pm4); 
        for i=1:NP
            if paint==1
                fit(i)=Obfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod,factor); 
            elseif paint==2
                fit(i)=SObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod); 
            elseif paint==3
                fit(i)=TObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod); 
            else
                fit(i)=FObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod);
            end
        end
        fitm=mean(fit);
        fitmax=max(fit);
        fitmin=min(fit);
        for i=1:NP
            if fit(i)<fitm
                CR(i)=CRL+(CRU-CRL)*(fitmax-fit(i))/(fitmax-fitmin);
            else
                CR(i)=CRL;
            end
        end
        CRM=[ ];
        for i=1:D
            CRM=[CRM,CR'];
        end
        mui=rand(NP,D)<CRM;
        mpo=mui < 0.5;
        ui=popold.*mpo+ui.*mui;
        for i=1:NP
             flag=DEtest(bound,ui(i,:));
             while flag==0
                 for j=1:D
                     ui(i,j)=pm5(i,j)+(2*rand-1)*(bound(j,2)...
                         -pm5(i,j));
                 end
                 flag=DEtest(bound,ui(i,:));
             end
        end
    end    
end
