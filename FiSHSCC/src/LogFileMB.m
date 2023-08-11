function LogFileMB(inputdata,outputfile,data,ScR,faultname,shearmodulusfrominput,straindropfrominput,M_dM_Lengths_forLOGoutput,out,n_sigma_trunc_magnitudes,yfc)
if ispc==1
  inputdata=  strrep(inputdata,'/','\');
end
fidout = fopen(strcat('./output_files/','LogFileMB', '.txt'), 'a');
fprintf(fidout, '*************************\n\n\n');
appending_pdf=fopen('./output_files/temporary_pdf_storage.txt');
% print the input file name followed by a blank line
fprintf(fidout,strcat('input file used:',blanks(1),inputdata,'\n'));
% print the output file name followed by a blank line
fprintf(fidout,strcat('output file created: ',blanks(1),outputfile,'\n'));
% print the time at which you did it followed by a blank line
t=clock;
fprintf(fidout,'executed at: ');
fprintf(fidout,'%d %d %d %d %d %3.2f\n\n',t);

% print the input data followed by a blank line
fprintf(fidout,strcat('input data used for calculation:\n'));
% print a title, followed by a blank line
fprintf(fidout, 'name ScaleRel length Dip Seismogenic_Thickness slipratemin slipratemax Mobs sdMobs occurrenceTime ShearMod(*10^10) StrainDrop(*10^-5) \n');
out1 = [data(:,2:9),shearmodulusfrominput,straindropfrominput];
for i=1:size(out1,1)
fprintf(fidout,'%s',faultname(i,:));
fprintf(fidout,'%1s',blanks(1));
fprintf(fidout,'%s',char(ScR(i,:)));
fprintf(fidout,'%1s',blanks(1));
fprintf(fidout,'%5.2f %5.2f %5.2f %5.2f %5.2f %3.1f %3.1f %d %5.3e %5.3e\n',out1(i,:));
end
fprintf(fidout,'\n\n');

%%% print the option used for caluclation
fprintf(fidout,strcat('number of sigma for trucated magnitudes function:',num2str(n_sigma_trunc_magnitudes),'\n'));
fprintf(fidout,strcat('year for calculating elasped time:',num2str(yfc),'\n\n'));


% print the output data followed by a blank line
fprintf(fidout,strcat('output data obtained:\n'));
% print a title, followed by a blank line
fprintf(fidout, 'name input-Lenght Lenght-from-AspectRatio MMo Mar Mscr1 Mscr2 Mobs dMMo dMar dMscr1 dMscr2 dMobs id Mmax sdMmax Tmean CV Telap Mo-rate\n');

out2=[M_dM_Lengths_forLOGoutput,out];
for i=1:size(out2,1)
fprintf(fidout,'%s',faultname(i));
fprintf(fidout,'%1s',blanks(1));
fprintf(fidout,'%5.2f %5.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %d %3.1f %3.1f %d %3.1f %d %5.3e\n',out2(i,:));
end

%% weigths(NaN if unweighted), x_range of magnitude and their unweighted-values
fprintf(fidout, '\n');
fprintf(fidout,'weigths(NaN if unweighted,0 if not used), x_range of magnitude and pdf_of_mag-values:\n\n');
tline = fgetl(appending_pdf);
fprintf(fidout, '%s\n', tline);
while ischar(tline)
        tline = fgetl(appending_pdf);
    
    if ischar(tline)
        fprintf(fidout, '%s\n', tline);
    end
    
end
fprintf(fidout, '*************************\n\n\n');
fclose(appending_pdf);
fclose(fidout);

