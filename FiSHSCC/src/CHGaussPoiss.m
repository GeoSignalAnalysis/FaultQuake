function [outputname]=CHGaussPoiss(c,d,outputfilename,faultname,mag,sdmag,Morate,id,nfault,w,Hpois,bin)

outputname=strcat(outputfilename,'_AR_ChGaussPoisson_rates', '.txt');
outputnameProbability=strcat(outputfilename,'_AR_ChGaussPoisson_Probability');   
  
% open two file for writing the outputs
fidout = fopen(strcat('./output_files/',outputname), 'w'); 
fidoutProb = fopen(strcat('./output_files/',outputnameProbability, '.txt'), 'w');
% print a title, followed by a blank line
fprintf(fidout, 'id Mmin bin rates name\n');
fprintf(fidoutProb, 'id Mmin window Probability name\n');





for i=1:nfault  % cycle for number of faults
  
magnitude_range=(mag(i)-sdmag(i)):bin:(mag(i)+sdmag(i));
M=10.^((c.*magnitude_range)+d);

pdf_mag = pdf('Normal',magnitude_range,mag(i),sdmag(i));
total_moment=sum(pdf_mag.*M);

ratio=(Morate(i))/total_moment;
balanced_pdf_moment=ratio*pdf_mag;

Mo_balanced(i,1)=sum(balanced_pdf_moment.*M); % for a  check

CHgaussRATES=balanced_pdf_moment;
cumCHgaussRATES=fliplr(cumsum(fliplr(CHgaussRATES)));

Mag_min=magnitude_range(1);
out_Rates=[id(i) Mag_min bin CHgaussRATES];
out_Prob = [id(i), Mag_min,w, Hpois(i)];

%%% PLOT figures and SAVE output files

fprintf(fidout,'%d, %3.1f, %3.1f,',out_Rates(1:3)); %id, minMag,bin
fprintf(fidout,'%1s',blanks(1));
for j=1:length(CHgaussRATES)
fprintf(fidout,'%5.4e',CHgaussRATES(j)); % rates per bin of the magnitude range 
fprintf(fidout,'%1s',blanks(1));
end

fprintf(fidout,',%1s',blanks(1));
fprintf(fidout,'%s\n',faultname(i,:));

fprintf(fidoutProb,'%d, %3.1f, %d, %5.3e,',out_Prob);
fprintf(fidoutProb,'%1s',blanks(1));
fprintf(fidoutProb,'%s\n',faultname(i,:));

figure(i)
semilogy((mag(i)-sdmag(i)):bin:(mag(i)+sdmag(i)),cumCHgaussRATES,'ok')
fault=faultname(i,:);
figname=strcat('./output_files/', outputfilename,'_AR_ChGaussPoisson_rates_',fault);

xlabel('magnitude');
ylabel('annual cumulative rates');
title(fault)
saveas(figure(i), figname,'epsc');

end