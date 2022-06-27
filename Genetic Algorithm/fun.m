function y=fun(x,reference,measure,N)
%% Fitness function
% reference   input: reference signal
% measure     input: measured signal
% N           input: the length of signal
cankao1=sigshift(reference,x(1),N);
y1=sigshift(cankao1,x(2),N);
y2=sigshift(cankao1,x(3),N);
y3=x(4)*y1+x(5)*y2;
y=sum((measure-y3).^2);
end





