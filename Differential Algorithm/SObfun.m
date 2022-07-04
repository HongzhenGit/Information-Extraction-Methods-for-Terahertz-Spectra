function y=SObfun(ff1,x,cc,m1,N,ml1,ml2,ml3,ml4,measure1,Freference1_obj,mod)
%% Fitness function for Double-layer samples
if mod==1
    e2=x(1)+(x(2)-x(1))./(1+1i*x(3)*ff1);
    e3=x(4)+(x(5)-x(4))./(1+1i*x(6)*ff1);
else
    e2=x(1)+x(2)^2./(x(3)^2-ff1.^2+1i*x(4)*ff1);
    e3=x(5)+x(6)^2./(x(7)^2-ff1.^2+1i*x(8)*ff1);
end
n1=1;
rn2=sqrt((abs(e2)+real(e2))/2);
ik2=sqrt((abs(e2)-real(e2))/2);
rn3=sqrt((abs(e3)+real(e3))/2);
ik3=sqrt((abs(e3)-real(e3))/2);
n2=rn2-1i*ik2;
n3=rn3-1i*ik3;
r12=(n2-n1)./(n1+n2);
r23=(n3-n2)./(n2+n3);                                                      
r34=1;
if mod==1
    Beta2=(2*pi*ff1.*n2*x(7)*10^(-2))/cc;  
    Beta3=(2*pi*ff1.*n3*x(8)*10^(-2))/cc;   
else
    Beta2=(2*pi*ff1.*n2*x(9)*10^(-2))/cc;  
    Beta3=(2*pi*ff1.*n3*x(10)*10^(-2))/cc; 
end
Hw3=(r23+r34.*exp(-2*1i*Beta3))./(1+r23.*r34.*exp(-2*1i*Beta3)); 
Hw2=(r12+Hw3.*exp(-2*1i*Beta2))./(1+r12.*Hw3.*exp(-2*1i*Beta2)); 
if mod==1
    SimuMeasure1=(Freference1_obj.*exp(2*1i*x(9)*10^(-2)*2*pi*ff1/cc)).*Hw2;
else
    SimuMeasure1=(Freference1_obj.*exp(2*1i*x(11)*10^(-2)*2*pi*ff1/cc)).*Hw2;
end
% The results of FFT is conjugated(Complex: a+bj, a-bj) 
for i=1:m1
    SimuMeasure2(i)=real(SimuMeasure1(m1-i+1))-1i*imag(SimuMeasure1(m1-i+1));
end
SimuMeasure=[zeros(1,ml1-1),SimuMeasure1,zeros(1,ml3-ml2-1),SimuMeasure2,zeros(1,N-ml4)];
TD_SimuMeasure=ifft(SimuMeasure);
y=sum((measure1-TD_SimuMeasure).^2);
end