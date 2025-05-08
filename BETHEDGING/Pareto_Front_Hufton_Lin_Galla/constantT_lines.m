%computes lines of constant T
A=2*max([kappa1,kappa2]);

compute_Tstar=1;
if compute_Tstar==1

%compute T*, the T corresponding to maximum growth rate.
pis_guess=[kappa1 kappa2]; %initial guess
        extra_param=[kappa1 kappa2 kA1  kB1 kA2 kB2 numpuntos_ave];
        disp('Maximizing growth rate')
        tic
       [optimal_pis,funcion_objetivo] = simulated_annealing('fun_lambda_for_minimization',pis_guess,extra_param)
%        [optimal_pis,funcion_objetivo] = simulated_annealing_var('fun_lambda_for_minimization',pis_guess,extra_param)
        toc
        disp('Computing Lambda for optimal pi''s')
	[V,lambda_max]=fun_variance_computation([kappa1 kappa2 optimal_pis(1) optimal_pis(2)],[ kA1  kB1 kA2 kB2],numpuntos_ave,numpuntos_var);

T_star=.5*sum(1./optimal_pis)
else
	T_star
end
%w=input('Pulsa una tecla')

pi1=[0.1:.02:3*A];
T=linspace(1/A,T_star,4)
%T=[0.5 1 2 T_star]
T=[0.5 .7 1.5 T_star]
Tlimits=T; %for fill_pareto
%T=1
%T=2/A
hold on
for t=T
	t
	pi2=1./(2*t-1./pi1)
	V=zeros(1,2*length(pi1));
	L=V;
	j=1
	for i=1:length(pi1)
		if pi2(i)>min(pi1) & pi2(i)<Inf
			pi2(i)
		[V(j),L(j)]=fun_variance_computation([kappa1 kappa2 pi1(i) pi2(i)],[ kA1  kB1 kA2 kB2],numpuntos_ave,numpuntos_var);
		j=j+1
		end
        end
   %     pi2=pi1;
%	pi1=1./(2*t-1./pi2)
%	for i=1:length(pi2)
%		if pi1(i)>0.3 & pi1(i)<Inf
%			pi1(i)
%		[V(j),L(j)]=fun_variance_computation([kappa1 kappa2 pi1(i) pi2(i)],[ kA1  kB1 kA2 kB2],numpuntos_ave,numpuntos_var);
%		j=j+1
%		end
 %       end
	std_dev=sqrt(V(1:j-1));
	L=L(1:j-1);
	plot(L,std_dev)
end


