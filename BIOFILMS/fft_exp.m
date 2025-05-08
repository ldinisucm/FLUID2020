function fft_num=fft_exp(tau,tc,w)
%Calcula la fft de una exponencial de tiempo tc, para cada tau (1-exp(-t/tc))*(exp(-t/tau))
%para las frecuencias w
for i=1:length(tau)
b=4;
M=2^10;
N=2*M;
deltaM=b/M;
tM=[1:deltaM:b];
h=exp(-tM/tau(i)).*(1-exp(-tM/tc));
h=[h zeros(1,N-M-1)];
fn=1/deltaM*[0:N/2-1]/N
