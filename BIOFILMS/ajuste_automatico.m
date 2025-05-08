%cargar los datos experimentales
tiempo=5;%12 24 48 72 96
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
Gp_exp=Gp_prune'
omega_exp_Gp=wp'
Gpp_exp=Gpp_prune'
omega_exp_Gpp=wpp'



%Guess inicial de los parámetros
eta=0.02;
x=1.01;
fc=3.5
Gp0=0.00005;
Gpp0_tix=0.000005;
Gpp0=2;
parametros=[x fc Gp0  eta Gpp0 Gpp0_tix];

tolerancia=1e-10
iter=500000
options=optimset('TolFun',tolerancia,'TolX',tolerancia,'MaxFunEvals',iter,'MaxIter',iter);
parmFit=fminsearch('fun_minimizar',parametros,options,omega_exp_Gp,omega_exp_Gpp,Gp_exp,Gpp_exp)

x=parmFit(1);
fc=parmFit(2);
wc=2*pi*fc;
Gp0=parmFit(3);
eta=parmFit(4);
Gpp0=parmFit(5);
Gpp0_tix=parmFit(6);

omega=[0.01 0.1 1:10 10:5:100 100:10:200];
X=omega./sinh(pi*omega./wc);
I0=besseli(0,X);
Gp_teor=Gp0*omega.^(x-1).*I0;

X=omega./sinh(pi*omega./wc);
I0=besseli(0,X);
Gpp_teor=Gpp0*(eta*omega+Gpp0_tix*omega.^(x-1).*I0);

loglog(omega,Gp_teor,omega,Gpp_teor)
