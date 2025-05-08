%Computes Pareto front from theoretical expressions for mean and variance

%Parameter definitions
set='a'


kappa_pi_set='fig2c'

switch kappa_pi_set
case {'fig2c'}
        kappa1 = 10; kappa2 = 10;
        %pi1=0.4;pi2=0.4;
	A=0.8
	beta0=10
	beta_final=6000
	n_iterations=2000;
pis_guess=[1 1]
case {'fig2b'}
  kappa1=1.0;kappa2=1.0;
  A=0.5
  beta0=10
  beta_final=1000
pis_guess=[kappa1 kappa2]; %initial guess
  %pi1=0.24;
  %pi2=pi1;
case {'fig2a'}
  kappa1=0.1;kappa2=0.1;
  %pi1=0.064;pi2=0.064;
case {'DL'}
kappa1 = 0.1; kappa2 = 0.33;
%pi1 = kappa1; pi2 = kappa2;
end

switch set
case {'a'}
%set a) Hufton
kA1 = 2; kB2 = -0.2;
kA2=-2; kB1=0.2;
case {'b'}
%set b) Hufton
kA1 = 0.5;
kB1=0.0001;
kA2=0.0001;
kB2=0.3250;
otherwise
        kA1=2;
        kB2=0.2;
        kA2=0;
        kB1=0;
end
%numpuntos_ave=1e6;
%numpuntos_var=1000;
numpuntos=4e4;

%kappa_set=[kappa1 kappa2];
%ks=[kA1 kB1 kA2 kB2];

Lambda=[];
Var=[];
PI0=[];
options=optimset();
%for alpha=.3:.05:1
alphas=[];
%for alpha=[0.99:-.01:0.95 0.9:-.05:.3]
for alpha=[.3:.025:.9 0.95:0.01:.99]
%for alpha=1
	alphas=[alphas alpha]
	extra_param=[kappa1 kappa2 kA1  kB1 kA2 kB2 alpha numpuntos];
	disp('Running minimization function')
	tic
  	[optimal_pis,funcion_objetivo] = simulated_annealing_var('fun_variance_computation_phi_for_Pareto',A,beta0,beta_final,n_iterations,pis_guess,extra_param)
	toc
	PI0=[PI0;optimal_pis];
	disp('Computing Lambda and VarLambda for optimal pi''s')
%	L=fun_lambda([kappa1 kappa2 optimal_pis(1) optimal_pis(2)],[ kA1  kB1 kA2 kB2],numpuntos_ave);
%	if alpha<1
%		V=(funcion_objetivo+alpha*L)/(1-alpha);
%	else
		[V,L]=fun_variance_computation_phi([kappa1 kappa2 optimal_pis(1) optimal_pis(2)],[ kA1  kB1 kA2 kB2],numpuntos);
%	end

	Lambda=[Lambda L];
	Var=[Var V];
	pis_guess=optimal_pis;

end
%sort befor plotting
[Lambda_sort,indices]=sort(Lambda);
Var_sort=Var(indices);
alpha_sort=alphas(indices);
PI0_1=PI0(indices,1);
PI0_2=PI0(indices,2);
plot(Lambda_sort,sqrt(Var_sort),'o-')
hold on
xlabel('Mean g.r.')
ylabel('std dev')
M=[alpha_sort' PI0_1 PI0_2 Lambda_sort' sqrt(Var_sort')]
save 'pareto_fig2c_phi.dat' -ascii  M
