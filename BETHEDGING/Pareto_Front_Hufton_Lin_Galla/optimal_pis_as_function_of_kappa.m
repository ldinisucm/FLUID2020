%Computes "Kelly's" like point

%Parameter definitions
set='a'

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
numpuntos_ave=2e6;
numpuntos_var=2000;
%numpuntos_ave=5e5;

%Annealing parameters
  beta0=1e3
  beta_final=1e7 %4000
n_iterations=100;
%pi1 = kappa1; pi2 = kappa2;
%kappa_set=[kappa1 kappa2];
%ks=[kA1 kB1 kA2 kB2];

Lambda=[];
Var=[];
PI0=[];
options=optimset();
%pis_guess=[0.3841    0.2777]
kappa1_vec=[1e-4 5e-4 1e-3 5e-3 0.01:.02:.09 0.1:.2:.9 1 1.1 1.2]
kappa1_vec=[1e-4 5e-4 1e-3 5e-3 0.01:.02:.09 0.1:.2:.9 .92 .96 1 1.02:.02:1.1]
%kappa1_vec=[.9 1 1.02:.02:1.1]
%kappa1_vec=[1e-4 5e-4 1e-3 5e-3]
kappa1_vec=[1e-3:3e-3:7e-3 0.01:.02:.09 0.1:0.05:.9]
        pis_guess=kappa1_vec(1)*[1 3.3];
num_kappas=length(kappa1_vec);
for i=1:num_kappas
	kappa1=kappa1_vec(i)
	kappa2=3.3*kappa1
	%A=kappa1/6
	A=1
	%pis_guess=[kappa1 kappa2]
	extra_param=[kappa1 kappa2 kA1  kB1 kA2 kB2 numpuntos_ave]; %for simulated annealing
%	extra_param=[kappa1 kappa2 kA1  kB1 kA2 kB2 numpuntos_ave]; %for fminsearch
	%disp('Running simulated annealing')
	disp('Running fminsearch')
	tic
  	%[optimal_pis,funcion_objetivo] = simulated_annealing_var('fun_lambda_for_minimization',A,beta0,beta_final,n_iterations,pis_guess,extra_param)
  	[optimal_pis,funcion_objetivo] = simulated_annealing_var_multiplicative('fun_lambda_for_minimization',A,beta0,beta_final,n_iterations,pis_guess,extra_param)
  	%[optimal_pis,funcion_objetivo] = simulated_annealing_var('fun_lambda_for_minimization_phi',A,beta0,beta_final,n_iterations,pis_guess,extra_param)
%	 [optimal_pis,funcion_objetivo] = fminsearch('fun_lambda_for_minimization_phi',pis_guess,options,extra_param)

	toc
	PI0=[PI0;optimal_pis];
	disp('Computing Lambda for optimal pi''s')
	%L=fun_lambda([kappa1 kappa2 optimal_pis(1) optimal_pis(2)],[ kA1  kB1 kA2 kB2],numpuntos_ave);
	[V,L]=fun_variance_computation([kappa1 kappa2 optimal_pis(1) optimal_pis(2)],[ kA1  kB1 kA2 kB2],numpuntos_ave,numpuntos_var);

	Lambda=[Lambda L];
	Var=[Var V];
        pis_guess=optimal_pis;
end
figure
%plot(kappa1_vec,PI0(:,1),kappa1_vec,PI0(:,2))
loglog(kappa1_vec,PI0(:,1),kappa1_vec,PI0(:,2))
ylabel('optimal pi')
xlabel('kappa1=kappa2/3.3')
legend('pi1','pi2')
figure
semilogx(kappa1_vec,Lambda,kappa1_vec,sqrt(Var))
xlabel('kappa1=kappa2/3.3')
legend('\Lambda','\sqrt{Var}')
