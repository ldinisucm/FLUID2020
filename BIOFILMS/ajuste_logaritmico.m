%cargar los datos experimentales
tiempo=4;%12 24 48 72 96
[omega_exp,Gp_exp,Gpp_exp]=datos_exp_fun_minimizar('Static',tiempo)
%Comprobamos si hay algún dato negativo y lo quitamos
Gp_prune=0;
wp=0;
Gpp_prune=0;
wpp=0;
jp=1;
jpp=1;
display('Quitando valores negativos')
for i=1:length(Gp_exp)
	if Gp_exp(i)>=0
		Gp_prune(jp)=Gp_exp(i);
		wp(jp)=omega_exp(i);
		jp=jp+1;
	end
end
for i=1:length(Gpp_exp)
        if Gpp_exp(i)>=0
		Gpp_prune(jpp)=Gpp_exp(i);
		wpp(jpp)=omega_exp(i);
		jpp=jpp+1;
	end
end
Gp_exp=log(Gp_prune)'
omega_exp_Gp=(wp)'
Gpp_exp=log(Gpp_prune)'
omega_exp_Gpp=(wpp)'



%Guess inicial de los parámetros
%Para 72h
eta=0.01;
x=1.1;
fc=6
Gp0=0.00002;
Gpp0_tix=0.000001;
Gpp0=1.2;
F=1;



%Bueno par 96h
%eta=0.01;
%x=1.1;
%fc=1.2
%Gp0=0.00002;
%Gpp0_tix=0.00001;
%Gpp0=1.2;
%F=4;
parametros=[x fc Gp0  eta Gpp0 Gpp0_tix F];

tolerancia=1e-8
iter=500000
options=optimset('TolFun',tolerancia,'TolX',tolerancia,'MaxFunEvals',iter,'MaxIter',iter);
options=optimset
parmFit=fminsearch('fun_minimizar_logaritmico',parametros,options,omega_exp_Gp,omega_exp_Gpp,Gp_exp,Gpp_exp)

x=parmFit(1);
fc=parmFit(2);
wc=2*pi*fc;
Gp0=parmFit(3);
eta=parmFit(4);
Gpp0=parmFit(5);
Gpp0_tix=parmFit(6);
F=parmFit(7);

omega=[0.01 0.1 1:10 10:5:100 100:10:200];
X=omega./sinh(pi*omega./wc);
I0=besseli(0,F*X);
Gp_teor=Gp0*omega.^(x-1).*I0;

X=omega./sinh(pi*omega./wc);
I0=besseli(0,F*X);
Gpp_teor=Gpp0*(eta*omega+Gpp0_tix*omega.^(x-1).*I0);

loglog(omega,Gp_teor,omega,Gpp_teor)
