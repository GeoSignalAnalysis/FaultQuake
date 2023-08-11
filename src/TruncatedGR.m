function [outputname]=TruncatedGR(c,d,outputfilename,faultname,mag,mt,Morate,id,nfault,bin,b)
global DATA

outputname=strcat(outputfilename,'_AR_TruncatedGR', '.txt');
  
fidout = fopen(strcat('./output_files/',outputname), 'w');
% print a title, followed by a blank line
fprintf(fidout, 'id Mmin bin rates name\n');





for i=1:nfault  % cycle for number of faults
  
magnitude_range=(mt(i):bin:(mag(i)));
M=10.^(c.*magnitude_range+d);
Beta=(2/3)*b(i);
Mt=10^(c*mt(i)+d);
Mxp=10^(c*(mag(i)+bin)+d); % added a bin to assume that the Max Magnitude given in input is reached.


TruncGR=((((Mt./M).^Beta) - ((Mt./Mxp).^Beta))./(1-((Mt./Mxp).^Beta))); 

Incremental=[fliplr(diff(fliplr(TruncGR))),TruncGR(end)];
Incremental_Morate=Incremental.*M;
Incremental_Morate_balanced=(Incremental_Morate.*Morate(i))./sum(Incremental_Morate);
cons_tassi_ind=(Incremental.*Incremental_Morate_balanced)./(Incremental_Morate);
cumulative_rates=fliplr(cumsum(fliplr(cons_tassi_ind)));
Mbalanced(i,1)=sum(cons_tassi_ind.*M);

out_Rates=[id(i) mt(i) bin cons_tassi_ind];

% DATA.(fault_name{i}).id=i;
DATA.(faultname{i}).rates=cons_tassi_ind;
DATA.(faultname{i}).bin=bin;





fprintf(fidout,'%d, %3.1f, %3.1f,',out_Rates(1:3));
fprintf(fidout,'%1s',blanks(1));
for j=1:length(cons_tassi_ind)
fprintf(fidout,'%5.4e',cons_tassi_ind(j));
fprintf(fidout,'%1s',blanks(1));
end

fprintf(fidout,',%1s',blanks(1));
fprintf(fidout,'%s\n',faultname{i});

figure(i)
semilogy(magnitude_range,cumulative_rates,'ok')
fault=faultname{i};
figname=strcat('./output_files/', outputfilename,'_AR_TruncatedGR_rates_',fault);

xlabel('magnitude');
ylabel('annual cumulative rates');
title(fault)
saveas(figure(i), figname,'epsc');
end
 export_faults_to_xml(DATA)
