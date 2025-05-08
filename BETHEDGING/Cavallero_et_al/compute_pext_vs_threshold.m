p=0.4
odds=[2.3 2]'
thres=[0:-0.01:-0.5];
num_races=17
pext=zeros(size(thres));
for i=1:length(thres)
	pext(i)=fun_prob_ext_computation_odds(p,odds,thres(i),num_races);
end
%figure
hold on
plot(thres,pext,'*-')
%xlabel('threshold')
%ylabel('pext')
