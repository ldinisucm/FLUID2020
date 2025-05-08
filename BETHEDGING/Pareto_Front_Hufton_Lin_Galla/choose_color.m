function color=choose_color(T,Tlimits)
%decides the color for a T
colores=['rygm'];
if length(Tlimits<=length(colores))
%colores=['rgbcmyk'];
else
	disp('Demasiados límites')
	cuantos=floor(Tlimits/length(colores))+1;
	for j=1:cuantos
		colores=[colores colores];
	end
end

color=colores(end);
for i=length(Tlimits):-1:1
	if T<Tlimits(i)
		color=colores(i);
	end
end
