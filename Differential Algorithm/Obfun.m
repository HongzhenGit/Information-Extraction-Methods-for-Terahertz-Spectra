function y=Obfun(ff1,x,cc,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod,factor)
%% Fitness function for single-lyaer samples
% ff1: input, the range of frequency spectrum
% x:   input, parameter vector
% cc:  input, the speed of light
% m1:  input, the lenth of ff1
% N:   input, the length of reference signal
% ml1,ml2,ml3,ml4: inputs, the position of spetrum(after cut)
% measure1:        input, measured signal in the time domain
% Freference1_obj: input, reference signal in the frequency domain
% mod:             input, the model indicator
% factor:          input, The number of mode in Debye or Lorentz model(1, 2)
if mod==1
    if factor==1
        e=x(1)+(x(2)-x(1))./(1+1i*x(3)*ff1);
    else
        e=x(1)+(x(2)-x(1))./(1+1i*x(3)*ff1)+(x(4)-x(2))./(1+1i*x(5)*ff1);
    end
else
    if factor==1
        e=x(1)+x(2)^2./(x(3)^2-ff1.^2+1i*x(4)*ff1);
    else
        e=x(1)+x(2)^2./(x(3)^2-ff1.^2+1i*x(4)*ff1)+x(5)^2./(x(6)^2-ff1.^2+1i*x(7)*ff1);
    end
end
n1=1;
ae=abs(e);
rn=sqrt((ae+real(e))/2);
ik=sqrt((ae-real(e))/2);
n2=rn-1i*ik;
r12=(n2-n1)./(n1+n2);
r23=1;                                                                    
if mod==1
    if factor==1
        Beta2=(2*pi*ff1.*n2*x(4)*10^(-2))/cc;                                      
    else
        Beta2=(2*pi*ff1.*n2*x(6)*10^(-2))/cc;
    end
else
    if factor==1
        Beta2=(2*pi*ff1.*n2*x(5)*10^(-2))/cc;                                      
    else
        Beta2=(2*pi*ff1.*n2*x(8)*10^(-2))/cc; 
    end
end
Hw=(r12+r23.*exp(-2*1i*Beta2))./(1+r12.*r23.*exp(-2*1i*Beta2));            
if mod==1
    if factor==1
        SimuMeasure1=(Freference1_obj.*exp(2*1i*x(5)*10^(-2)*2*pi*ff1/cc)).*Hw;
    else
        SimuMeasure1=(Freference1_obj.*exp(2*1i*x(7)*10^(-2)*2*pi*ff1/cc)).*Hw;
    end
else
    if factor==1
        SimuMeasure1=(Freference1_obj.*exp(2*1i*x(6)*10^(-2)*2*pi*ff1/cc)).*Hw;
    else
        SimuMeasure1=(Freference1_obj.*exp(2*1i*x(9)*10^(-2)*2*pi*ff1/cc)).*Hw;
    end
end
% The results of FFT is conjugated(Complex: a+bj, a-bj) 
for i=1:m1
    SimuMeasure2(i)=real(SimuMeasure1(m1-i+1))-1i*imag(SimuMeasure1(m1-i+1));
end
SimuMeasure=[zeros(1,ml1-1),SimuMeasure1,zeros(1,ml3-ml2-1),SimuMeasure2,...
    zeros(1,N-ml4)];
TD_SimuMeasure=ifft(SimuMeasure);
y=sum((measure1-TD_SimuMeasure).^2);
end