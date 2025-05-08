%cargar los datos experimentales
tiempo=5;%12 24 48 72 96
datos_exp_fun('Static',tiempo)

%omega=[1:10 10:5:100 100:10:200 200:50:300];
omega=[0.01 0.1 1:10 10:5:100 100:10:200];


x=1.01
wc=1
tau0=1
a=0;
atc=00;
a_dft=0;
aheavi=0;
b=0.2
wc=0
N=100000;
Gprima_num=0;
dft_num=0;
omega=10.^[-2:.01:2];
G0=4e-2;
eta0=.01
%Gpp0=0.1

tau1=25
tau2=0.0
G1=0
for i=1:length(omega)
  w=omega(i);
  %u=[tau0*w:Du:10000*tau0*w];
  tau=linspace(tau0,10000*tau0,N);
  %tau2=linspace(tau0,100*tau0,N);
  u=linspace(tau0*w,10000*tau0*w,N);
  Du=(max(u)-min(u))/N;
  Dtau=(max(tau)-min(tau))/N;
  Gast=(j*w*tau-tau2*tau*w^2)./(1+(1+G1)*j*w*tau+j*w*tau2-w^2*tau*(tau1+tau2));
  int_tau=Dtau*sum(tau.^(-x).*real(Gast));
  gpp(i)=G0*(Dtau*sum(tau.^(-x).*imag(Gast)))+eta0*w;

int1=Du*sum(u.^(2-x)./(1+u.^2));
%int2=Du*sum(u.^(1-x)./(1+u.^2));
Z=1/G0;
a(i)=int1*w^(x-1)/Z;
atc(i)=int_tau/Z;
%aheavi(i)=int_heavi/Z;
end
loglog(omega,a,omega,atc)
loglog(omega,gpp)
