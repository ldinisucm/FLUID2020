function [optimal_pis,funcion_objetivo] = simulated_annealing(nombre_funcion_para_minimizar,pis_guess,extra_param)
%minimizes nombre_funcion.
%nombre_funcion has to admit pis_guess as input, and some other parameters

tol=1e-4;
diff=tol+1;
xcurr=pis_guess;
n_iterations=100
beta0=10;
beta_final=1000;
J0=feval(nombre_funcion_para_minimizar,xcurr,extra_param);
J0=real(J0)
tamx=size(xcurr);
for n=0:n_iterations
	n
% beta=beta0*(1-n/n_iterations)+n/n_iterations*beta_final	
 logbeta=log(beta0)*(1-n/n_iterations)+n/n_iterations*log(beta_final)
 beta=exp(logbeta) 
 A=1/sqrt(beta);
 A=0.5
 xnew=xcurr+(rand(tamx)-0.5)*A
 xcurr
 %check for valid xnew (all elements positive)
 while any(xnew<0)
         disp('Negativo!')
         xnew=xcurr+(rand(tamx)-0.5)*A
 end
 Jnew=feval(nombre_funcion_para_minimizar,xnew,extra_param);
 Jnew=real(Jnew)
 J0
 resta=Jnew-J0
 if Jnew<J0
	 %accept
	 J0=Jnew;
	 xcurr=xnew;
	 disp('--------------Aceptado por menor------------------')
 else
	 random_number=rand();
         exp_factor=exp(-beta*(Jnew-J0))
	 if random_number < exp(-beta*(Jnew-J0))
	       %accept
	 disp('Aceptado por exp')
              J0=Jnew;
              xcurr=xnew;
	 else
		 disp('---no se acepta---')
            end
 end
end
optimal_pis=xcurr
funcion_objetivo=J0


end

