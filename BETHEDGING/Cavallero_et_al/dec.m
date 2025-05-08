function index=dec(vector)
%Takes a binary vector and returns the proper index which is the decimal+1. Units are the last cipher, and, the last but one is 2^1, etc..
N=length(vector);
power=[2.^(N-1:-1:0)];
%power=[2.^(0:1:N-1)];
index=sum(power.*vector)+1;

