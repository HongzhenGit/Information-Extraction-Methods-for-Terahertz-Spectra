function y=sigshift(x,m,N)
%% Shift operation
% y=sigshift(x,m,N)
%y£ºthe output sequence after shift
%x£ºinput sequence
%m£ºthe number of shifts
%N£ºthe lentght of input sequence
if length(x)>N
    error('N must be >= the length of x')    
end
x=[x zeros(1,N-length(x))];
n=[0:1:N-1];
n=mod(n-m,N);
y=x(n+1);