
function [outputname]=ClassicalGR(c,d,outputfilename,faultname,mag,mt,Morate,id,nfault,bin,b)

outputname=strcat(outputfilename,'_AR_ClassicalGR', '.txt');
    
fidout = fopen(strcat('./output_files/',outputname), 'w');
% print a title, followed by a blank line
fprintf(fidout, 'id Mmin bin rates name\n');



for i=1:nfault  % cycle for number of faults
a=1; % starting value for a  
magnitude_range=(mt(i):bin:(mag(i)));
M=10.^(c.*magnitude_range+d);
  
Ncum=10.^(a-b(i).*magnitude_range);
Incremental=[fliplr(diff(fliplr(Ncum))),Ncum(end)];
Incremental_Morate=Incremental.*M;
Incremental_Morate_balanced=(Incremental_Morate.*Morate(i))./sum(Incremental_Morate);
cons_tassi_ind=(Incremental.*Incremental_Morate_balanced)./(Incremental_Morate);
cumulative_rates=fliplr(cumsum(fliplr(cons_tassi_ind)));
Mbalanced(i,1)=sum(cons_tassi_ind.*M);

out_Rates=[id(i) mt(i) bin cons_tassi_ind];

fprintf(fidout,'%d, %3.1f, %3.1f,',out_Rates(1:3));
fprintf(fidout,'%1s',blanks(1));
for j=1:length(cons_tassi_ind)
fprintf(fidout,'%5.4e',cons_tassi_ind(j));
fprintf(fidout,'%1s',blanks(1));
end

fprintf(fidout,',%1s',blanks(1));
fprintf(fidout,'%s\n',faultname(i,:));

figure(i)
semilogy(magnitude_range,cumulative_rates,'ok')
fault=faultname(i,:);
figname=strcat('./output_files/', outputfilename,'_AR_ClassicalGR_rates_',fault);

xlabel('magnitude');
ylabel('annual cumulative rates');
title(fault)
saveas(figure(i), figname,'epsc');
end