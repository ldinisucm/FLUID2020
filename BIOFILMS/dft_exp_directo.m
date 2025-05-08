function dft_num=dft_exp_directo(tau,tc,w)
%Calcula la transformada de fourier para una w concreta haciendo la integral directamente, lo que no es muy bueno, pero es fácil
tf=1e3*tau/w;
N=5e7;
Dt=tf/N;
t=[0:Dt:tf];
dft_num=0;
for i=1:length(tau)
  h=exp(-(t/tau(i))).*(1-exp(-t/tc));
  dft_num(i)=Dt*sum(h.*exp(-j*w*t));
end
