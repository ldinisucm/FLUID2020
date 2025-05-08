function dft_num=dft_exp_directo_vector(tau,tc,w)
%Calcula la transformada de fourier para una w concreta haciendo la integral directamente, lo que no es muy bueno, pero es fácil
tf=1000*tau/w;
N=1e5;
Dt=tf/N;
dft_num=0;
wtc=.4;
for k=1:length(tf)
  t=[0:Dt:tf(k)];
  %h=exp(-(t/tau(k))).*(1-exp(-t/tc));
  h=exp(-(t/tau(k))).*tanh((t-tc)/wtc);
  dft_num(k)=Dt(k)*sum(h.*exp(-j*w*t));
end
