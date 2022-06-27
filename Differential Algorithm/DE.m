%% DE Algorithm
clc
clear all
close all
mod=1;               % The model for Reflection Index: Debye(1) or Lorentz(2)                           
factor=1;            % The number of mode in Debye or Lorentz model(1, 2)
paint=1;             % The number of layers of sample(Up to 4)
strategy=1;          % The strategy for Corss/Mutation, please check DECrossMutation.m for more details
itermax=200;         % The maximum number of iterations
%% Initializing
% Parameter Bound
Lb1=1;                                                                      
Hb1=5;
Lb4=100;                                                                     
Hb4=400;
Lb5=200;                                                                      
Hb5=500;
if mod==1
    Lb2=1;                                                                     
    Hb2=10;
    Lb3=0;                                                                      
    Hb3=15;
    if paint==1
        if factor==1
            D=5;
            NP=70;
            bound=[Lb1,Hb1;Lb2,Hb2;Lb3,Hb3;Lb4,Hb4;Lb5,Hb5];
        else
            D=7;
            NP=70;
            bound=[Lb1,Hb1;Lb2,Hb2;Lb3,Hb3;Lb2,Hb2;Lb3,Hb3;Lb4,Hb4;Lb5,Hb5];
        end
    elseif paint==2
        D=9;
        NP=200;
        bound=[Lb1,Hb1;Lb2,Hb2;Lb3,Hb3;Lb1,Hb1;Lb2,Hb2;Lb3,Hb3;Lb4,Hb4;...
               Lb4,Hb4;Lb5,Hb5];
    elseif paint==3
        D=13;
        NP=130;
        bound=[Lb1,Hb1;Lb2,Hb2;Lb3,Hb3;Lb1,Hb1;Lb2,Hb2;Lb3,Hb3;Lb1,Hb1;...
               Lb2,Hb2;Lb3,Hb3;Lb4,Hb4;Lb4,Hb4;Lb4,Hb4;Lb5,Hb5];
    else
        D=17;
        NP=170;
        bound=[Lb1,Hb1;Lb2,Hb2;Lb3,Hb3;Lb1,Hb1;Lb2,Hb2;Lb3,Hb3;Lb1,Hb1;...
               Lb2,Hb2;Lb3,Hb3;Lb1,Hb1;Lb2,Hb2;Lb3,Hb3;Lb4,Hb4;Lb4,Hb4;...
               Lb4,Hb4;Lb4,Hb4;Lb5,Hb5];
    end
else
    Lb2=0;                                                                     
    Hb2=200;
    Lb3=0;                                                                      
    Hb3=500;
    if paint==1
        if factor==1
            D=6;
            NP=60;
            bound=[Lb1,Hb1;Lb2,Hb2;Lb2,Hb2;Lb3,Hb3;Lb4,Hb4;Lb5,Hb5];
        else
            D=9;
            NP=90;
            bound=[Lb1,Hb1;Lb2,Hb2;Lb2,Hb2;Lb3,Hb3;Lb2,Hb2;Lb2,Hb2;...
                   Lb3,Hb3;Lb4,Hb4;Lb5,Hb5];
        end
    elseif paint==2
        D=11;
        NP=200;
        bound=[Lb1,Hb1;Lb2,Hb2;Lb2,Hb2;Lb3,Hb3;Lb1,Hb1;Lb2,Hb2;Lb2,Hb2;...
               Lb3,Hb3;Lb4,Hb4;Lb4,Hb4;Lb5,Hb5];
    elseif paint==3
        D=16;
        NP=160;
        bound=[Lb1,Hb1;Lb2,Hb2;Lb2,Hb2;Lb3,Hb3;Lb1,Hb1;Lb2,Hb2;Lb2,Hb2;...
               Lb3,Hb3;Lb1,Hb1;Lb2,Hb2;Lb2,Hb2;Lb3,Hb3;Lb4,Hb4;Lb4,Hb4;...
               Lb4,Hb4;Lb5,Hb5];
    else
        D=21;
        NP=210;
        bound=[Lb1,Hb1;Lb2,Hb2;Lb2,Hb2;Lb3,Hb3;Lb1,Hb1;Lb2,Hb2;Lb2,Hb2;...
               Lb3,Hb3;Lb1,Hb1;Lb2,Hb2;Lb2,Hb2;Lb3,Hb3; Lb1,Hb1;Lb2,Hb2;...
               Lb2,Hb2;Lb3,Hb3;Lb4,Hb4;Lb4,Hb4;Lb4,Hb4;Lb4,Hb4;Lb5,Hb5];
    end
end                                                                        
                                                                  
%% Input Checking
if (NP < 5)
   fprintf(1,'Error! NP should be >= 5\n');
end
if (itermax < 0)
   fprintf(1,'Error! itermax should be > 0\n');
end
% Construct the initial population
pop=DECode(NP,D,bound);  

%% Matrix to store information
popold=zeros(size(pop)); % Population Matrix during interation
val=zeros(1,NP);         % Matrix for fitness value of each individual

%% Read and process the signals
MeaData=textread('20190524干燥过程10：20.txt');                                        
RefData=textread('20190524干燥过程参考.txt');                                        
ts=RefData(:,1)';
dt=ts(2)-ts(1);
N=length(ts);
% Get the frequency axis
n=0:N-1;
f=n/(N*dt);                                                                
c=3;  
measure=10*MeaData(:,2);
refrence=10*RefData(:,2);
%% EEMD+Wavelet: Remove non-flat baseline and noises
% Wavelet filtering for reference signal
wname  = 'sym4';               % Wavelet for analysis.
sorh   = 's';                  % Type of thresholding.
level = 1:1:20;                % Level for wavelet decomposition.
nb_Int = 1:1:6;                % Number of intervals for thresholding.
SNR_refence = zeros(20,6);
for i = 1:20
    for j = 1:6
        referenceDen = cmddenoise(refrence,wname,level(i),sorh,nb_Int(j));
        SNR_refence(i,j) = SNR_singlech(referenceDen',refrence);
    end
end
[row,col] = find(SNR_refence==min(min(SNR_refence)));
referenceDen_wave = cmddenoise(refrence,wname,level(row),sorh,nb_Int(col));
referenceDen_wave = referenceDen_wave';

% EEMD Signal Re-construction for reference signal
emd_reference = eemd(referenceDen_wave,0.2,150);
reference_nobase = zeros(length(referenceDen_wave),1);                                       
for i=3:7
    reference_nobase = reference_nobase + emd_reference(:,i);
end
SNR1_reference = zeros(20,6);
for i = 1:20
    for j = 1:6
        referenceDen1 = cmddenoise(reference_nobase,wname,level(i),sorh,nb_Int(j));
        SNR1_reference(i,j) = SNR_singlech(referenceDen1',reference_nobase);
    end
end
[row1,col1] = find(SNR1_reference==min(min(SNR1_reference)));
referenceDen_emd = cmddenoise(reference_nobase,wname,level(row1),sorh,nb_Int(col1));

% Wavelet Filtering for measured signal
SNR_measure = zeros(20,6);
for i = 1:20
    for j = 1:6
        measureDen = cmddenoise(measure,wname,level(i),sorh,nb_Int(j));
        SNR_measure(i,j) = SNR_singlech(measureDen',measure);
    end
end
[row,col] = find(SNR_measure==min(min(SNR_measure)));
measureDen_wave = cmddenoise(measure,wname,level(row),sorh,nb_Int(col));
measureDen_wave = measureDen_wave';

% EEMD signal Re-construction for measured signal
emd_measure = eemd(measureDen_wave,0.2,150);
measure_nobase = zeros(length(measureDen_wave),1);                                       
for i=3:7
    measure_nobase = measure_nobase + emd_measure(:,i);
end
SNR1_measure = zeros(20,6);
for i = 1:20
    for j = 1:6
        measureDen1 = cmddenoise(measure_nobase,wname,level(i),sorh,nb_Int(j));
        SNR1_measure(i,j) = SNR_singlech(measureDen1',measure_nobase);
    end
end
[row1,col1] = find(SNR1_measure==min(min(SNR1_measure)));
measureDen_emd = cmddenoise(measure_nobase,wname,level(row1),sorh,nb_Int(col1));
refrence = referenceDen_emd;
measure = measureDen_emd;

%% Fourier Tranformation and Spectrum Cut
Fmeasure1=fft(measure);
Frefrence1=fft(refrence);                                                    
ff1=f-0.1;                                                                 
ff2=f-2.5;
[mv1,ml1]=min(abs(ff1));
[mv2,ml2]=min(abs(ff2));
f1=f(ml1:ml2);
m1=length(f1);
Freference1_obj=Frefrence1(ml1:ml2);
Fmeasure1_obj=Fmeasure1(ml1:ml2);
pos1=find(abs(Fmeasure1)==abs(Fmeasure1_obj(1)));
pos2=find(abs(Fmeasure1)==abs(Fmeasure1_obj(m1)));
ml3=pos2(2);
ml4=pos1(2);
Frefrence1=[zeros(1,ml1-1),Frefrence1(ml1:ml2),zeros(1,ml3-ml2-1),...
    Frefrence1(ml3:ml4),zeros(1,N-ml4)]; 
Fmeasure1=[zeros(1,ml1-1),Fmeasure1(ml1:ml2),zeros(1,ml3-ml2-1),...
    Fmeasure1(ml3:ml4),zeros(1,N-ml4)]; 
refrence=ifft(Frefrence1);
measure=ifft(Fmeasure1);

%% Fitness-Initial Evaluation
for i=1:NP
    if paint==1
        val(i)=Obfun(f1,pop(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure,Freference1_obj,mod,factor);
    elseif paint==2
        val(i)=SObfun(f1,pop(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure,Freference1_obj,mod);
    elseif paint==3
        val(i)=TObfun(f1,pop(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure,Freference1_obj,mod); 
    else
        val(i)=FObfun(f1,pop(i,:),c,m1,N,ml1,ml2,ml3,ml4,measure,Freference1_obj,mod); 
    end
end
[bestval,bestindex]=min(val);
bestmemit=pop(bestindex,:);
bestmem=bestmemit;
nfeval=NP;

%% Initialize random matrics
pm1=zeros(NP,D);                                                           % initialize population matrix 1
pm2=zeros(NP,D);                                                           % initialize population matrix 2
pm3=zeros(NP,D);                                                           % initialize population matrix 3
pm4=zeros(NP,D);                                                           % initialize population matrix 4
pm5=zeros(NP,D);                                                           % initialize population matrix 5
bm=zeros(NP,D);                                                            % initialize bestmember  matrix
ui=zeros(NP,D);                                                            % intermediate population of perturbed vectors
mui=zeros(NP,D);                                                           % mask for intermediate population
mpo=zeros(NP,D);                                                           % mask for old population
rot=(0:1:NP-1);                                                            % rotating index array
rt=zeros(NP);                                                              % another rotating index array
a1=zeros(NP);                                                              % index array1
a2=zeros(NP);                                                              % index array2
a3=zeros(NP);                                                              % index array3
a4=zeros(NP);                                                              % index array4
a5=zeros(NP);                                                              % index array5
ind=zeros(4);
iter = 1;
bb(1)=bestval;

%% Evolution start
while ((iter<itermax)&&(bestval>1.e-6))   
  popold=pop;                                                              % save the old population

  ui=DECrossMutation(strategy,rot,NP,D,popold,bestmemit,bound,f1,c,m1,N,...
      ml1,ml2,ml3,ml4,measure,Freference1_obj,mod,paint,factor);           % Cross and Mutation

  
  for i=1:NP
      if paint==1
         tempval=Obfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,...
             measure,Freference1_obj,mod,factor);                          % check cost of competitor
      elseif paint==2
         tempval=SObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,...
             measure,Freference1_obj,mod);                                 
      elseif paint==3
        tempval=TObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,...
             measure,Freference1_obj,mod);                                 
      else          
        tempval=FObfun(f1,ui(i,:),c,m1,N,ml1,ml2,ml3,ml4,...
             measure,Freference1_obj,mod);                                 
      end
    nfeval=nfeval + 1;
    if (tempval<=val(i))                                                   % Only keep the one with better fitness
       pop(i,:)=ui(i,:);                                                   
       val(i)=tempval;                                                                                                                              
       if (tempval<bestval)                                                
          bestval=tempval;                                                 % new best value
          bestmem=ui(i,:);                                                 % new best parameter vector ever
       end
    end
  end                                                                      
  bestmemit=bestmem;                                                       % freeze the best member of this iteration for the coming 
  iter = iter + 1;                                                         % iteration. This is needed for some of the strategies.
  bb(iter)=bestval;
  iter
end
%% Results
if paint==1
    %% For Single Layer
    if mod==1
        if factor==1
            model=bestmem(1:3);
            d1=bestmem(4);
            d=bestmem(5);
        else
            model=bestmem(1:5);
            d1=bestmem(6);
            d=bestmem(7);
        end
    else
        if factor==1
            model=bestmem(1:4);
            d1=bestmem(5);
            d=bestmem(6);
        else
            model=bestmem(1:7);
            d1=bestmem(8);
            d=bestmem(9);
        end
    end
    disp '**************Debye/Lorentz Model Parameters****************'
    model
    disp '*************************Thickness**************************'
    d1
    disp '********************Reference Position**********************'
    d
    min(bb)
    % Calculate the time domain signal
    if mod==1
        if factor==1
            e=bestmem(1)+(bestmem(2)-bestmem(1))./(1+1i*bestmem(3)*f1);
        else
            e=bestmem(1)+(bestmem(2)-bestmem(1))./(1+1i*bestmem(3)*f1)+(bestmem(4)-bestmem(2))./(1+1i*bestmem(5)*f1);
        end
    else
        if factor==1
            e=bestmem(1)+bestmem(2)^2./(bestmem(3)^2-f1.^2+1i*bestmem(4)*f1);
        else
            e=bestmem(1)+bestmem(2)^2./(bestmem(3)^2-f1.^2+1i*bestmem(4)*f1)+bestmem(5)^2./(bestmem(6)^2-f1.^2+1i*bestmem(7)*f1);
        end
    end
    n1=1;
    ae=abs(e);
    rn=sqrt((ae+real(e))/2);
    ik=sqrt((ae-real(e))/2);
    n2=rn-1i*ik;
    r12=(n2-n1)./(n1+n2);
    % Fresnel Number
    r23=1;                                                                     
    if mod==1
        % Phase Shift   
        if factor==1
            Beta2=(2*pi*f1.*n2*bestmem(4)*10^(-2))/c;                                  
        else
            Beta2=(2*pi*f1.*n2*bestmem(6)*10^(-2))/c;
        end
    else
        if factor==1
            Beta2=(2*pi*f1.*n2*bestmem(5)*10^(-2))/c;                                
        else
            Beta2=(2*pi*f1.*n2*bestmem(8)*10^(-2))/c; 
        end
    end
    % Rouard's method for single-layer medium
    Hw=(r12+r23.*exp(-2*1i*Beta2))./(1+r12.*r23.*exp(-2*1i*Beta2));            
    if mod==1
        if factor==1
            SimMeasure1=(Freference1_obj.*exp(2*1i*bestmem(5)*10^(-2)*2*pi*f1/c)).*Hw;
        else
            SimMeasure1=(Freference1_obj.*exp(2*1i*bestmem(7)*10^(-2)*2*pi*f1/c)).*Hw;
        end
    else
        if factor==1
            SimMeasure1=(Freference1_obj.*exp(2*1i*bestmem(6)*10^(-2)*2*pi*f1/c)).*Hw;
        else
            SimMeasure1=(Freference1_obj.*exp(2*1i*bestmem(9)*10^(-2)*2*pi*f1/c)).*Hw;
        end
    end
    for i=1:m1
        SimMeasure2(i)=real(SimMeasure1(m1-i+1))-1i*imag(SimMeasure1(m1-i+1));
    end
    SimMeasure=[zeros(1,ml1-1),SimMeasure1,zeros(1,ml3-ml2-1),SimMeasure2,...
        zeros(1,N-ml4)];
    TDsimMeasure=ifft(SimMeasure);
    % Dielectric Constant and Reflection Index
    figure
    plot(f1,real(e),'r-','LineWidth',2)
    hold on
    plot(f1,imag(e),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Dielectric Constant')
    figure
    plot(f1,real(n2),'r-','LineWidth',2)
    hold on
    plot(f1,imag(n2),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Reflection Index')
elseif paint==2
    %% For Double Layer
    if mod==1
        lorentz2=bestmem(1:3);
        lorentz3=bestmem(4:6);
        d2=bestmem(7);
        d3=bestmem(8);
        d=bestmem(9);
    else
        lorentz2=bestmem(1:4);
        lorentz3=bestmem(5:8);
        d2=bestmem(9);
        d3=bestmem(10);
        d=bestmem(11);
    end
    disp '************Debye/Lorentz model parameters(Surface)**************'
    lorentz2
    disp '**********************Thickness(Surface)*************************'
    d2
    disp '************Debye/Lorentz model parameters(Bottom)***************'
    lorentz3
    disp '***********************Thickness(Bottom)*************************'
    d3
    disp '**********************Reference Position*************************'
    d
    disp '*********************Final Fitted Residual***********************'
    min(bb)
    % Calculate Time Domain Signal
    if mod==1
        e2=bestmem(1)+(bestmem(2)-bestmem(1))./(1+1i*bestmem(3)*f1);
        e3=bestmem(4)+(bestmem(5)-bestmem(4))./(1+1i*bestmem(6)*f1);
    else
        e2=bestmem(1)+bestmem(2)^2./(bestmem(3)^2-f1.^2+1i*bestmem(4)*f1);
        e3=bestmem(5)+bestmem(6)^2./(bestmem(7)^2-f1.^2+1i*bestmem(8)*f1);
    end
    n1=1;
    rn2=sqrt((abs(e2)+real(e2))/2);
    ik2=sqrt((abs(e2)-real(e2))/2);
    rn3=sqrt((abs(e3)+real(e3))/2);
    ik3=sqrt((abs(e3)-real(e3))/2);
    n2=rn2-1i*ik2;
    n3=rn3-1i*ik3;
    % Fresnel Number
    r12=(n2-n1)./(n1+n2);
    r23=(n3-n2)./(n2+n3);                                                      
    r34=1;
    if mod==1
        Beta2=(2*pi*f1.*n2*bestmem(7)*10^(-2))/c;  
        Beta3=(2*pi*f1.*n3*bestmem(8)*10^(-2))/c;   
    else
        Beta2=(2*pi*f1.*n2*bestmem(9)*10^(-2))/c;  
        Beta3=(2*pi*f1.*n3*bestmem(10)*10^(-2))/c;   
    end
    % Rouard's method for double-layer medium
    Hw3=(r23+r34.*exp(-2*1i*Beta3))./(1+r23.*r34.*exp(-2*1i*Beta3)); 
    Hw2=(r12+Hw3.*exp(-2*1i*Beta2))./(1+r12.*Hw3.*exp(-2*1i*Beta2)); 
    if mod==1
        SimMeasure1=(Freference1_obj.*exp(2*1i*bestmem(9)*10^(-2)*2*pi*f1/c)).*Hw2;
    else
        SimMeasure1=(Freference1_obj.*exp(2*1i*bestmem(11)*10^(-2)*2*pi*f1/c)).*Hw2;
    end
    for i=1:m1
        SimMeasure2(i)=real(SimMeasure1(m1-i+1))-1i*imag(SimMeasure1(m1-i+1));
    end
    SimMeasure=[zeros(1,ml1-1),SimMeasure1,zeros(1,ml3-ml2-1),SimMeasure2,zeros(1,N-ml4)];
    TDsimMeasure=ifft(SimMeasure);
    % Dielectric Constant and Reflection Index
    figure
    plot(f1,real(e2),'--','LineWidth',2)
    hold on
    plot(f1,imag(e2),'k-','LineWidth',2)
    legend('Real','Imaginary')
    title('Di-Constant(Surface)')
    figure
    plot(f1,real(e3),'--','LineWidth',2)
    hold on
    plot(f1,imag(e3),'k-','LineWidth',2)
    legend('Real','Imaginary')
    title('Di-Constant(Bottom)')
    figure
    plot(f1,real(n2),'--','LineWidth',2)
    hold on
    plot(f1,imag(n2),'k-','LineWidth',2)
    legend('Real','Imaginary')
    title('Ref-Index(Surface)')
    figure
    plot(f1,real(n3),'--','LineWidth',2)
    hold on
    plot(f1,imag(n3),'k-','LineWidth',2)
    legend('Real','Imaginary')
    title('Ref-Index(Bottom)')
elseif paint==3
    %% For tripple layer
    if mod==1
        lorentz2=bestmem(1:3);
        lorentz3=bestmem(4:6);
        lorentz4=bestmem(7:9);
        d2=bestmem(10);
        d3=bestmem(11);
        d4=bestmem(12);
        d=bestmem(13);
    else
        lorentz2=bestmem(1:4);
        lorentz3=bestmem(5:8);
        lorentz4=bestmem(9:12);
        d2=bestmem(13);
        d3=bestmem(14);
        d4=bestmem(15);
        d=bestmem(16);
    end
    disp '************Debye/Lorentz Model Parameters(Surface)**************'
    lorentz2
    disp '***********************Thickness(Surface)************************'
    d2
    disp '************Debye/Lorentz Model Parameters(Medium)***************'
    lorentz3
    disp '**********************Thickness(Medium)**************************'
    d3
    disp '************Debye/Lorentz Model Parameters(Bottom)***************'
    lorentz4
    disp '**********************Thickness(Bottom)**************************'
    d4
    disp '**********************Reference Position*************************'
    d
    disp '*********************Final Fitted Residual***********************'
    min(bb)
    % Calculate Time Domain Signal
    if mod==1
        e2=bestmem(1)+(bestmem(2)-bestmem(1))./(1+1i*bestmem(3)*f1);
        e3=bestmem(4)+(bestmem(5)-bestmem(4))./(1+1i*bestmem(6)*f1);
        e4=bestmem(7)+(bestmem(8)-bestmem(7))./(1+1i*bestmem(9)*f1);
    else
        e2=bestmem(1)+bestmem(2)^2./(bestmem(3)^2-f1.^2+1i*bestmem(4)*f1);
        e3=bestmem(5)+bestmem(6)^2./(bestmem(7)^2-f1.^2+1i*bestmem(8)*f1);
        e4=bestmem(9)+bestmem(10)^2./(bestmem(11)^2-f1.^2+1i*bestmem(12)*f1);
    end
    n1=1;
    rn2=sqrt((abs(e2)+real(e2))/2);
    ik2=sqrt((abs(e2)-real(e2))/2);
    rn3=sqrt((abs(e3)+real(e3))/2);
    ik3=sqrt((abs(e3)-real(e3))/2);
    rn4=sqrt((abs(e4)+real(e4))/2);
    ik4=sqrt((abs(e4)-real(e4))/2);
    n2=rn2-1i*ik2;
    n3=rn3-1i*ik3;
    n4=rn4-1i*ik4;
    % Fresnel Number
    r12=(n2-n1)./(n1+n2);
    r23=(n3-n2)./(n2+n3);  
    r34=(n4-n3)./(n3+n4);
    r45=1;
    if mod==1
        Beta2=(2*pi*f1.*n2*bestmem(10)*10^(-2))/c;  
        Beta3=(2*pi*f1.*n3*bestmem(11)*10^(-2))/c;   
        Beta4=(2*pi*f1.*n4*bestmem(12)*10^(-2))/c; 
    else
        Beta2=(2*pi*f1.*n2*bestmem(13)*10^(-2))/c;  
        Beta3=(2*pi*f1.*n3*bestmem(14)*10^(-2))/c;   
        Beta4=(2*pi*f1.*n4*bestmem(15)*10^(-2))/c; 
    end
    % Rouard's method for tripple-layer medium
    Hw4=(r34+r45.*exp(-2*1i*Beta4))./(1+r34.*r45.*exp(-2*1i*Beta4)); 
    Hw3=(r23+Hw4.*exp(-2*1i*Beta3))./(1+r23.*Hw4.*exp(-2*1i*Beta3)); 
    Hw2=(r12+Hw3.*exp(-2*1i*Beta2))./(1+r12.*Hw3.*exp(-2*1i*Beta2)); 
    if mod==1
        SimMeasure1=(Freference1_obj.*exp(2*1i*bestmem(13)*10^(-2)*2*pi*f1/c)).*Hw2;
    else
        SimMeasure1=(Freference1_obj.*exp(2*1i*bestmem(16)*10^(-2)*2*pi*f1/c)).*Hw2;
    end
    for i=1:m1
        SimMeasure2(i)=real(SimMeasure1(m1-i+1))-1i*imag(SimMeasure1(m1-i+1));
    end
    SimMeasure=[zeros(1,ml1-1),SimMeasure1,zeros(1,ml3-ml2-1),SimMeasure2,zeros(1,N-ml4)];
    TDsimMeasure=ifft(SimMeasure);
    % Dielectric Constant and Reflection Index
    figure
    plot(f1,real(e2),'r-','LineWidth',2)
    hold on
    plot(f1,imag(e2),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Di-Constant(Surface)')
    figure
    plot(f1,real(e3),'r-','LineWidth',2)
    hold on
    plot(f1,imag(e3),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Di-Constant(Medium)')
    figure
    plot(f1,real(e4),'r-','LineWidth',2)
    hold on
    plot(f1,imag(e4),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Di-Constant(Bottom)')
    figure
    plot(f1,real(n2),'r-','LineWidth',2)
    hold on
    plot(f1,imag(n2),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Ref-Index(Surface)')
    figure
    plot(f1,real(n3),'r-','LineWidth',2)
    hold on
    plot(f1,imag(n3),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Ref-Index(Medium)')
    figure
    plot(f1,real(n4),'r-','LineWidth',2)
    hold on
    plot(f1,imag(n4),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Ref-Index(Bottom)')
else
    %% For four layer 
    if mod==1
        model2=bestmem(1:3);
        model3=bestmem(4:6);
        model4=bestmem(7:9);
        model5=bestmem(10:12);
        d2=bestmem(13);
        d3=bestmem(14);
        d4=bestmem(15);
        d5=bestmem(16);
        d=bestmem(17);
    else
        model2=bestmem(1:4);
        model3=bestmem(5:8);
        model4=bestmem(9:12);
        model5=bestmem(13:16);
        d2=bestmem(17);
        d3=bestmem(18);
        d4=bestmem(19);
        d5=bestmem(20);
        d=bestmem(21);
    end
    disp '*************Debye/Lorentz Model Parameters(First)***************'
    model2
    disp '************************Thickness(First)*************************'
    d2
    disp '************Debye/Lorentz Model Parameters(Second)***************'
    model3
    disp '************************Thickness(Second)************************'
    d3
    disp '*************Debye/Lorentz Model Parameters(Third)***************'
    model4
    disp '************************Thickness(Third)*************************'
    d4
    disp '************Debye/Lorentz Model Parameters(Fourth)***************'
    model5
    disp '***********************Thickness(Fourth)*************************'
    d5
    disp '**********************Reference Position*************************'
    d
    disp '*********************Final Fitted Residual***********************'
    min(bb)
    % Calculate Time Domain Signal
    if mod==1
        e2=bestmem(1)+(bestmem(2)-bestmem(1))./(1+1i*bestmem(3)*f1);
        e3=bestmem(4)+(bestmem(5)-bestmem(4))./(1+1i*bestmem(6)*f1);
        e4=bestmem(7)+(bestmem(8)-bestmem(7))./(1+1i*bestmem(9)*f1);
        e5=bestmem(10)+(bestmem(11)-bestmem(10))./(1+1i*bestmem(12)*f1);
    else
        e2=bestmem(1)+bestmem(2)^2./(bestmem(3)^2-f1.^2+1i*bestmem(4)*f1);
        e3=bestmem(5)+bestmem(6)^2./(bestmem(7)^2-f1.^2+1i*bestmem(8)*f1);
        e4=bestmem(9)+bestmem(10)^2./(bestmem(11)^2-f1.^2+1i*bestmem(12)*f1);
        e5=bestmem(13)+bestmem(14)^2./(bestmem(15)^2-f1.^2+1i*bestmem(16)*f1);
    end
    n1=1;
    rn2=sqrt((abs(e2)+real(e2))/2);
    ik2=sqrt((abs(e2)-real(e2))/2);
    rn3=sqrt((abs(e3)+real(e3))/2);
    ik3=sqrt((abs(e3)-real(e3))/2);
    rn4=sqrt((abs(e4)+real(e4))/2);
    ik4=sqrt((abs(e4)-real(e4))/2);
    rn5=sqrt((abs(e5)+real(e5))/2);
    ik5=sqrt((abs(e5)-real(e5))/2);
    n2=rn2-1i*ik2;
    n3=rn3-1i*ik3;
    n4=rn4-1i*ik4;
    n5=rn5-1i*ik5;
    % Fresnel Number
    r12=(n2-n1)./(n1+n2);
    r23=(n3-n2)./(n2+n3);                                                      
    r34=(n4-n3)./(n3+n4);
    r45=(n5-n4)./(n4+n5);
    r56=1;
    if mod==1
        Beta2=(2*pi*f1.*n2*bestmem(13)*10^(-2))/c;  
        Beta3=(2*pi*f1.*n3*bestmem(14)*10^(-2))/c;   
        Beta4=(2*pi*f1.*n4*bestmem(15)*10^(-2))/c; 
        Beta5=(2*pi*f1.*n5*bestmem(16)*10^(-2))/c; 
    else
        Beta2=(2*pi*f1.*n2*bestmem(17)*10^(-2))/c;  
        Beta3=(2*pi*f1.*n3*bestmem(18)*10^(-2))/c;   
        Beta4=(2*pi*f1.*n4*bestmem(19)*10^(-2))/c; 
        Beta5=(2*pi*f1.*n5*bestmem(20)*10^(-2))/c; 
    end
    % Rouard's method for four-layer medium
    Hw5=(r45+r56.*exp(-2*1i*Beta5))./(1+r45.*r56.*exp(-2*1i*Beta5)); 
    Hw4=(r34+Hw5.*exp(-2*1i*Beta4))./(1+r34.*Hw5.*exp(-2*1i*Beta4)); 
    Hw3=(r23+Hw4.*exp(-2*1i*Beta3))./(1+r23.*Hw4.*exp(-2*1i*Beta3)); 
    Hw2=(r12+Hw3.*exp(-2*1i*Beta2))./(1+r12.*Hw3.*exp(-2*1i*Beta2)); 
    if mod==1
        SimMeasure1=(Freference1_obj.*exp(2*1i*bestmem(17)*10^(-2)*2*pi*f1/c)).*Hw2;
    else
        SimMeasure1=(Freference1_obj.*exp(2*1i*bestmem(21)*10^(-2)*2*pi*f1/c)).*Hw2;
    end
    for i=1:m1
        SimMeasure2(i)=real(SimMeasure1(m1-i+1))-1i*imag(SimMeasure1(m1-i+1));
    end
    SimMeasure=[zeros(1,ml1-1),SimMeasure1,zeros(1,ml3-ml2-1),SimMeasure2,zeros(1,N-ml4)];
    TDsimMeasure=ifft(SimMeasure);
    % Dielectric Constant and Reflection Index
    figure
    plot(f1,real(e2),'r-','LineWidth',2)
    hold on
    plot(f1,imag(e2),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Di-Constant(First)')
    figure
    plot(f1,real(e3),'r-','LineWidth',2)
    hold on
    plot(f1,imag(e3),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Di-Constant(Second)')
    figure
    plot(f1,real(e4),'r-','LineWidth',2)
    hold on
    plot(f1,imag(e4),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Di-Constant(Third)')
    figure
    plot(f1,real(e5),'r-','LineWidth',2)
    hold on
    plot(f1,imag(e5),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Di-Constant(Fourth)')
    figure
    plot(f1,real(n2),'r-','LineWidth',2)
    hold on
    plot(f1,imag(n2),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Ref-Index(First)')
    figure
    plot(f1,real(n3),'r-','LineWidth',2)
    hold on
    plot(f1,imag(n3),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Ref-Index(Second)')
    figure
    plot(f1,real(n4),'r-','LineWidth',2)
    hold on
    plot(f1,imag(n4),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Ref-Index(Third)')
    figure
    plot(f1,real(n5),'r-','LineWidth',2)
    hold on
    plot(f1,imag(n5),'b-','LineWidth',2)
    legend('Real','Imaginary')
    title('Ref-Index(Fourth)')
end
%% Fitness Curve
figure
plot(bb,'LineWidth',2)
title(['Fitness Curve' 'Iteration Stopped＝' num2str(itermax)]);
xlabel('Number of iterations');ylabel('Fitness');
grid on
%% Some Plotting Commands
figure
plot(f1,abs(SimMeasure(ml1:ml2)),'r-','LineWidth',2)
hold on
plot(f1,abs(Fmeasure1_obj),'k-','LineWidth',2)
legend('Fitted','Actual')
title('Amplitude(Frequency Domain)')
figure
plot(f,unwrap(angle(SimMeasure)),'r-','LineWidth',2)
hold on
plot(f,unwrap(angle(Fmeasure1)),'k-','LineWidth',2)
legend('Fitted','Actual')
title('Phase(Frequency Domain)')
figure
plot(ts,TDsimMeasure,'r-','LineWidth',2)
hold on
plot(ts,measure,'k-','LineWidth',2)
legend('Fitted','Actual')
title('Time Domain Signal')