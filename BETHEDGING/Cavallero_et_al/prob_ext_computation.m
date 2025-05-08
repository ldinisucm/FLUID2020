logC0=0;
tray=[0;1]
num_pasos=14
p=.6;
P=[1-p;p]

odds=[2.3 2]';
b=P;
delta_plus=log(b(1)*odds(1));
delta_minus=log(b(2)*odds(2));


delta_plus=0.6;
delta_minus=-0.6;
threshold=-1;
if logC0+delta_minus<threshold
	P(1)=0;
end
if logC0+delta_plus<threshold
	P(2)=0;
end
P

for i=2:num_pasos
	aux=[tray;tray];
	col_nueva=[zeros(length(tray),1);ones(length(tray),1)];
	tray=[col_nueva aux]
	[f,c]=size(tray);
	Pold=P
	P=zeros(f,1);
	for j=1:f
		tray(j,:)	
		logC=logC0+sum(tray(j,:))*delta_plus-(c-sum(tray(j,:)))*delta_minus
 		if logC>threshold
			index=dec(tray(j,:))
			previous_index=dec(tray(j,1:i-1))
			P(index)=(p*tray(j,i)+(1-p)*(1-tray(j,i)))*Pold(previous_index);
			P(index)
		end
	end
	P
end
Pext=1-sum(P)
