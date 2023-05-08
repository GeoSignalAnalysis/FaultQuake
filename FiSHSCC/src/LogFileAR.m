function LogFileAR(inputdata,outputname,data,faultname,Fault_behaviour,MFD,window, binstep, grinputoptions,ParametersGR)
if ispc==1
  inputdata=  strrep(inputdata,'/','\');
end
fidout = fopen(strcat('./output_files/','LogFileAR', '.txt'), 'a');
fprintf(fidout, '*************************\n\n\n');
% print the input file name followed by a blank line
fprintf(fidout,strcat('input file used:',blanks(1),inputdata,'\n'));
% print the output file name followed by a blank line
fprintf(fidout,strcat('output file created: ',blanks(1),outputname,'\n'));
% print the time at which you did it followed by a blank line
t=clock;
fprintf(fidout,'executed at: ');
fprintf(fidout,'%d %d %d %d %d %3.2f\n\n',t);

% print the input data followed by a blank line
fprintf(fidout,strcat('input data used for calculation:\n'));
% print a title, followed by a blank line
% in the case of models 1,2,4,5,7,8 the input file format is 
%id Mchar sdMchar Tmean alfa Telap Mo_rate(N/m2) name
if Fault_behaviour == 1 |   Fault_behaviour == 2 | Fault_behaviour == 4 | Fault_behaviour == 5 | Fault_behaviour == 7 |  Fault_behaviour == 8
  fprintf(fidout, 'id Mchar sdMchar Tmean alfa Telap Mo_rate(N/m2) name\n');
  for i=1:size(data,1)
fprintf(fidout,'%d %5.2f %5.2f %d %5.2f %d %5.3e',data(i,:));
fprintf(fidout,'%1s',blanks(1));
fprintf(fidout,'%s\n',faultname(i,:));
end
% in the case of models 3,6 the input file format is 
%id Mchar sdMchar Tmean alfa Telap Mo_rate(N/m2) Probability name
elseif Fault_behaviour == 3 |   Fault_behaviour == 6
fprintf(fidout, 'id Mchar sdMchar Tmean alfa Telap Mo_rate(N/m2) Probability name\n');
  for i=1:size(data,1)
fprintf(fidout,'%d %5.2f %5.2f %d %5.2f %d %5.3e %5.3e',data(i,:));
fprintf(fidout,'%1s',blanks(1));
fprintf(fidout,'%s\n',faultname(i,:));
  end
end

fprintf(fidout,'\n\n');

%%% print the option used for caluclation
fprintf(fidout,strcat('Magnitude-Frequency-Distribution:',MFD,'\n'));
if ~isnan(grinputoptions)
fprintf(fidout,'GR options (id, b-value, Mt):\n');
fprintf(fidout,'%d %3.2f %3.2f\n',ParametersGR');
end
fprintf(fidout,strcat('forecast period:',num2str(window),'\n'));
fprintf(fidout,strcat('binstep:',num2str(binstep),'\n\n'));

fprintf(fidout, '*************************\n\n\n');
fclose(fidout);

