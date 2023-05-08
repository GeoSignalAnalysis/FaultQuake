function LogFileRP(name,outputfilename,data,seedrnd,number_of_simulations)
if ispc==1
  name=  strrep(name,'/','\');
end
fidout = fopen(strcat('./output_files/','LogFileRP', '.txt'), 'a');
fprintf(fidout, '*************************\n\n\n');
% print the input file name followed by a blank line
fprintf(fidout,strcat('input file used:',blanks(1),name,'\n'));
% print the output file name followed by a blank line
fprintf(fidout,strcat('output file created: ',blanks(1),outputfilename,'\n'));
% print the time at which you did it followed by a blank line
t=clock;
fprintf(fidout,'executed at: ');
fprintf(fidout,'%d %d %d %d %d %3.2f\n\n',t);

% print the input data followed by a blank line
fprintf(fidout,strcat('input data used for calculation:\n'));
% print a title, followed by a blank line
fprintf(fidout, 'upper-limit lower-limit\n');
out1 = [data];
fprintf(fidout,'%d %d\n',out1');

fprintf(fidout,'\n\n');

%%% print the option used for caluclation
fprintf(fidout,strcat('number of simulations:',num2str(number_of_simulations),'\n'));
if ~isnan(seedrnd)
fprintf(fidout,strcat('you have specified the seed:',num2str(seedrnd),'\n\n'));
else
fprintf(fidout,strcat('you have not specified the seed\n\n'));
end
fprintf(fidout, '*************************\n\n\n');
fclose(fidout);

