function [outputname]=singleProb(outputfilename,mag,id,nfault,faultname,w,ProbabilityOfOccurrence)


outputname=(strcat(outputfilename,'_AR_SingleProb_rates', '.txt'));
% open a file for writing the output
fidout = fopen(strcat('./output_files/',outputname), 'w');
% print a title, followed by a blank line
fprintf(fidout, 'id Mchar rate name\n');




%%% calculate a fictious Tmean following Pace et al., 2006
Tfict=(-1*w)./log(1-ProbabilityOfOccurrence);


CH_RATES(:,1)=1./Tfict;
CH_MAGNITUDE(:,1)=mag;

out_Rates = [id, CH_MAGNITUDE,CH_RATES];





%%% PLOT figures and SAVE output files
for i =1:nfault
fprintf(fidout,'%d, %3.1f, %5.3e,',out_Rates(i,:));
fprintf(fidout,'%1s',blanks(1));
fprintf(fidout,'%s\n',faultname(i,:));

figure(i)
semilogy(out_Rates(i,2),out_Rates(i,3),'ok')
fault=faultname(i,:);
figname=strcat('./output_files/', outputfilename,'_AR_SingleProb_rates_',fault);

xlabel('magnitude');
ylabel('annual rate');
title(fault)
saveas(figure(i), figname,'epsc');

end