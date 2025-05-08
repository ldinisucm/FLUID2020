%cargar los datos experimentales
tiempo=4;%12 24 48 72 96
datos_exp_fun('Static',tiempo)

%omega=[1:10 10:5:100 100:10:200 200:50:300];
omega=[0.01 0.1 1:10 10:5:100 100:10:200];


x=1.01
wc=1
tau0=10
a=0;
atc=00;
a_dft=0;
aheavi=0;
tc=50;
N=1000;
G0=1;
Gprima_num=0;
dft_num=0;
omega=10.^[-2:.32:2];
for i=1:length(omega)
  w=omega(i);
  %u=[tau0*w:Du:10000*tau0*w];
  tau=linspace(tau0,100*tau0,N);
  u=linspace(tau0*w,100*tau0*w,N);
  Du=(max(u)-min(u))/N;
  Dtau=(max(tau)-min(tau))/N;
  %int_tau=Dtau*sum((w^2*tau.^(2-x))./(1+w^2*tau.^2));
  int_tau=Dtau*sum(heaviside(tau-tc).*(w^2*tau.^(2-x))./(1+w^2*tau.^2));
 % int_heavi=Dtau*sum(exp(-tc./tau).*(w^2*tau.^(2-x).*cos(w*tc)+w*tau.^(1-x)*sin(w*tc))./(1+w^2*tau.^2));
int1=Du*sum(u.^(2-x)./(1+u.^2));
%int2=Du*sum(u.^(1-x)./(1+u.^2));
Z=1/G0;
a(i)=int1*w^(x-1)/Z;
atc(i)=int_tau/Z;
%aheavi(i)=int_heavi/Z;
end
loglog(omega,a,omega,atc)
