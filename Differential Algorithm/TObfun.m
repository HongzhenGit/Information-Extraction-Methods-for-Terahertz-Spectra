function y=TObfun(ff1,x,cc,m1,N,ml1,ml2,ml3,ml4,celiang1,fcankao1obj,mod)
%% Fitness function for Tripple-layer samples
if mod==1
    e2=x(1)+(x(2)-x(1))./(1+1i*x(3)*ff1);
    e3=x(4)+(x(5)-x(4))./(1+1i*x(6)*ff1);
    e4=x(7)+(x(8)-x(7))./(1+1i*x(9)*ff1);
else
    e2=x(1)+x(2)^2./(x(3)^2-ff1.^2+1i*x(4)*ff1);
    e3=x(5)+x(6)^2./(x(7)^2-ff1.^2+1i*x(8)*ff1);
    e4=x(9)+x(10)^2./(x(11)^2-ff1.^2+1i*x(12)*ff1);
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
r12=(n2-n1)./(n1+n2);
r23=(n3-n2)./(n2+n3);                                                      
r34=(n4-n3)./(n3+n4);
r45=1;
if mod==1
    Beta2=(2*pi*ff1.*n2*x(10)*10^(-2))/cc;  
    Beta3=(2*pi*ff1.*n3*x(11)*10^(-2))/cc;   
    Beta4=(2*pi*ff1.*n4*x(12)*10^(-2))/cc; 
else
    Beta2=(2*pi*ff1.*n2*x(13)*10^(-2))/cc;  
    Beta3=(2*pi*ff1.*n3*x(14)*10^(-2))/cc;   
    Beta4=(2*pi*ff1.*n4*x(15)*10^(-2))/cc; 
end
Hw4=(r34+r45.*exp(-2*1i*Beta4))./(1+r34.*r45.*exp(-2*1i*Beta4)); 
Hw3=(r23+Hw4.*exp(-2*1i*Beta3))./(1+r23.*Hw4.*exp(-2*1i*Beta3)); 
Hw2=(r12+Hw3.*exp(-2*1i*Beta2))./(1+r12.*Hw3.*exp(-2*1i*Beta2)); 
if mod==1
    Simceliang1=(fcankao1obj.*exp(2*1i*x(13)*10^(-2)*2*pi*ff1/cc)).*Hw2;
else
    Simceliang1=(fcankao1obj.*exp(2*1i*x(16)*10^(-2)*2*pi*ff1/cc)).*Hw2;
end
for i=1:m1
    Simceliang2(i)=real(Simceliang1(m1-i+1))-1i*imag(Simceliang1(m1-i+1));
end
Simceliang=[zeros(1,ml1-1),Simceliang1,zeros(1,ml3-ml2-1),Simceliang2,zeros(1,N-ml4)];
TDsimceliang=ifft(Simceliang);
y=sum((celiang1-TDsimceliang).^2);
end