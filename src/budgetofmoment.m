%%% BUDGET OF MOMENT
%%% CODE TO COMPUTE Tmean and CV from geometry and slip rates data



global imfile
global outputname
global Weighted
global WeightInputFile

global StrainDropFromInput
global StrainDrop
global ShearModulusFromInput
global ShearModulus
global n_sigma
global year_for_calculations
global MB_output


ShearModulusFromInputFile=ShearModulusFromInput;
StrainDropFromInputFile=StrainDropFromInput;
if isempty(ShearModulusFromInputFile)
    ShearModulusFromInputFile=0;
end
if isempty(StrainDropFromInputFile)
    StrainDropFromInputFile=0;
end
%%% WARNING MESSAGES
warning off backtrace
warning off verbose
if isempty(ShearModulus)  & ShearModulusFromInputFile~=1;
  warning('Shear modulus: You are using default value 3*10^10 Pa')
 ShearModulus=3;
end
if isempty(StrainDrop) & StrainDropFromInputFile~=1;
  warning('Strain Drop: You are using default value 3*10^-5')
 StrainDrop=3;
end
if isempty(n_sigma) | isnan(n_sigma)
    warning('you are using untruncated gaussian distributions for magnitudes')
    n_sigma_trunc_magnitudes=0;
else
    n_sigma_trunc_magnitudes=n_sigma;
end
if isempty(year_for_calculations)
  warning('Year: You are using current year for computing elapsed time')
 year_for_calculations=clock;yfc=year_for_calculations(1);
else
    yfc=year_for_calculations;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open a file for writing the output
outputfile=strcat('./output_files/',outputname,'_MB', '.txt');
fidout = fopen(outputfile, 'w');

% print a title on the output file, followed by a blank line
fprintf(fidout, 'id Mmax sdMmax Tmean CV Telap Mo-rate name\n');

temporary_pdf_storage=fopen('./output_files/temporary_pdf_storage.txt', 'w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('RUNNING...\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% DEFINITION OF COEFFICIENTS
%%%% coefficients by Hanks & Kanamori, 1985
d=9.1; c=1.5;
dMMO=0.3; %standard deviation of the magnitude calculated by the seismic moment definition
%%%%%%%%%%%%%%%%%%
%%% IF PARAMETERS ARE FROM OPTIONS, not from input file
if ShearModulusFromInputFile~=1;
mu_opt=ShearModulus*10^10;
end
if StrainDropFromInputFile~=1;
straindrop_opt=StrainDrop*10^-5;
end
%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% LOAD INPUT FILE
%%% 4 cases: 
%%% 1) ShearModulus and StrainDrop not given (10 coloumns);
%%% 2) ShearModulus in the input file but StrainDrop not given (11 coloumns);
%%% 3) StrainDrop in the input file but ShearModulus not given (11 coloumns);
%%% 4) ShearModulus and StrainDrop in the input file (12 coloumns);

%%% CASE 1
if (ShearModulusFromInputFile~=1) & (StrainDropFromInputFile~=1)
inputdata=imfile;
fid = fopen(inputdata);
head1=fscanf(fid, '%s %s %s %s %s %s %s %s %s %s %s\n',11);
F=textscan(fid,'%s %s %f %f %f %f %f %f %f %f %f');
fault_name=F{1};
faultname=char(fault_name);
nfault=size(faultname(:,:),1);
ScaleRel=F{2};
data=[(1:nfault)',F{3:10}];
SCC = F{11};
shearmodulusfrominput(1:size(data,1),1)=mu_opt;
straindropfrominput(1:size(data,1),1)=straindrop_opt;
end
%%% CASE 2
if (ShearModulusFromInputFile==1) & (StrainDropFromInputFile~=1)
inputdata=imfile;
fid = fopen(inputdata);
head1=fscanf(fid, '%s %s %s %s %s %s %s %s %s %s %s\n',12);
F=textscan(fid,'%s %s %f %f %f %f %f %f %f %f %f %f');
fault_name=F{1};
faultname=char(fault_name);
nfault=size(faultname(:,:),1);
ScaleRel=F{2};
data=[(1:nfault)',F{3:10}];
shearmodulusfrominput=[F{11}]*10^10;
straindropfrominput(1:size(data,1),1)=straindrop_opt;
end
%%% CASE 3
if (ShearModulusFromInputFile~=1) & (StrainDropFromInputFile==1)
inputdata=imfile;
fid = fopen(inputdata);
head1=fscanf(fid, '%s %s %s %s %s %s %s %s %s %s %s\n',11);
F=textscan(fid,'%s %s %f %f %f %f %f %f %f %f %f');
fault_name=F{1};
faultname=char(fault_name);
nfault=size(faultname(:,:),1);
ScaleRel=F{2};
data=[(1:nfault)',F{3:10}];
shearmodulusfrominput(1:size(data,1),1)=mu_opt;
straindropfrominput=[F{11}]*10^-5;
end

%%% CASE 4
if (ShearModulusFromInputFile==1) & (StrainDropFromInputFile==1)
inputdata=imfile;
fid = fopen(inputdata);
head1=fscanf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s\n',12);
F=textscan(fid,'%s %s %f %f %f %f %f %f %f %f %f %f');
fault_name=F{1};
faultname=char(fault_name);
nfault=size(faultname(:,:),1);
ScaleRel=F{2};
data=[(1:nfault)',F{3:10}];
shearmodulusfrominput=[F{11}]*10^10;
straindropfrominput=[F{12}]*10^-5;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:nfault; % count the number of faults

%% read input parameters
code=data(i,1); 
ScR=ScaleRel(i,:);
Length=data(i,2);
Dip=data(i,3); 
Seismogenic_thickness=data(i,4); 
Slipmin=data(i,5); 
Slipmax=data(i,6); 
mag=data(i,7); 
sdmag=data(i,8); 
Telap=yfc-data(i,9);
fault=faultname(i,:);

mu=shearmodulusfrominput(i);
straindrop=straindropfrominput(i);
%%%%%%%%%%%%%%%% weigths
%%% LOAD WEIGHTS for PDFs,if given as a file
 if ~isempty(Weighted)
     inputWeight=WeightInputFile;
     weights_fault=load(inputWeight);
 end
%%%%%%%%%%%%%%%%%%

% convert from km to m, and from mm/yr to m/yr
Length=Length*1000;
Width=(Seismogenic_thickness*1000)/sind(Dip);
V=(Slipmin+Slipmax)/2000;
dV=V-(Slipmin/1000);
% return the appropriate coefficients of the used scale-relationships
[coeff,ARtable]=kin2coeff(ScR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% MAIN BLOCK 1  calculate Mmax with a standard deviation    %%%%

% return the appropriate magnitudes from the used scale-relationships 
[MRLD MRA dMRLD dMRA MAR ar_coeff LAR legends_Mw]=coeff2mag(ScR,coeff,Length,Width,ARtable,mu,straindrop);

% calculate the moment magnitude
MMO= (1/c)* (log10(straindrop * mu* Length^2 *Width) -d);

% create a vector of magnitudes M and a vector of sdandard deviations dM
M=[MMO; MAR; MRLD; MRA; mag];
M=round(M*100)/100; %round to second decimal
dM=[dMMO; ar_coeff(3); dMRLD; dMRA; sdmag];
dM=round(dM*100)/100; %round to second decimal

M_dM_Lengths_forLOGoutput(i,:)=[Length/1000,LAR/1000, M',dM'];

M(isnan(M))=[]; % if a magnitude is equal to NaN then it is not used
dM(isnan(dM))=[]; % if a magnitude is equal to NaN then it is not used


% calculate Mmax
clear xmag, clear y, clear sumy
clear muhat, clear sigmahat, clear sdm

%define a vector of common range of magnitude with a bin 0.01 for computing pdf from vectors M and dM
x_range_of_mag=[];
if n_sigma_trunc_magnitudes==0
x_range_of_mag=floor(min(M-dM)):0.01:ceil(max(M+dM)); 
else
x_range_of_mag=floor(min(M-n_sigma_trunc_magnitudes*dM)):0.01:ceil(max(M+n_sigma_trunc_magnitudes*dM));
end

% for each M calculate a probability density function using a normal distribution
for k=1:length(M)
    pdf_magnitudes(k,:)=normpdf(x_range_of_mag,M(k),dM(k));
end



% for the Scale-Relationships derived magnitude & observed Magnitude, if desired, calculate the
% truncated distribution given the number of sd
if n_sigma_trunc_magnitudes>0
[dist_trunc_mag]=truncGaussDist(pdf_magnitudes,x_range_of_mag,M,dM,mag,n_sigma_trunc_magnitudes);   
pdf_magnitudes(1:(size(dist_trunc_mag,1)),:)=dist_trunc_mag;
end

% normalize each pdf to its maximum value, so each pdf weigths as 1.
pdf_magnitudes=pdf_magnitudes./repmat((max(pdf_magnitudes,[],2)),1,size(pdf_magnitudes,2));

%% check that the  size of weights and magnitudes is the same
%% if obs magnitude does not exist then recompute weigths to have their sum=1
%% if LAR >= Length its weigth is 0 
if ~isempty(Weighted)
  if LAR>=Length
    weights_fault(2)=0;
  end
weights=weights_fault(1:length(M))/sum(weights_fault(1:length(M)));
pdf_magnitudes=pdf_magnitudes.*repmat((weights),1,size(pdf_magnitudes,2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% weighted pdf are stored and saved in the LogFileMB with their weights
out_pdf=[];out_weights=[];
if isempty(Weighted)
    out_weights(1:length(M),1)=NaN;
else
    out_weights(1:length(M),1)=weights;
end

if LAR>=Length
    out_weights(2)=0;
end

out_x_range=strcat(['weights/magnitude',blanks(1),num2str(x_range_of_mag)]);

out_pdf=[out_weights,pdf_magnitudes];

fprintf(temporary_pdf_storage,'%s\n',faultname(i,:));
fprintf(temporary_pdf_storage,out_x_range);
fprintf(temporary_pdf_storage,'\n');
for tps=1:size(out_pdf,1)
    fprintf(temporary_pdf_storage,strcat(num2str(out_pdf(tps,:)),'\n'));
end
fprintf(temporary_pdf_storage,'\n\n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%


% calcuate the summed distribution
% if LAR is greater or equal than Length then LAR is not used 
if LAR>=Length
    pdf_magnitudes(2,:)=[];
end
summed_pdf_magnitudes=sum(pdf_magnitudes);

% calculate mean and standard deviation of the summed distribution by
% fitting it with a normal distribution
[Mmax,sigma_Mmax] = normfit(x_range_of_mag,[],[],summed_pdf_magnitudes); 
Mmax=round(Mmax*10)/10; % round to the first decimal
sigma_Mmax=round(sigma_Mmax*10)/10; % round to the first decimal

% save the Mmax and its standard deviation 
output_Mmax_sigmaMmax(i,:)=[Mmax,sigma_Mmax];

% from Mmax and sigma_Mmax calculate a probability density function using a normal distribution
% actually, this is not in the graph
gauss_fit=normpdf(x_range_of_mag,Mmax,sigma_Mmax);
if ~isempty(Weighted)
gauss_fit=gauss_fit/max(gauss_fit);
end
%%%% END OF MAIN BLOCK 1                                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN BLOCK 2 plots                                        %%%%

%%%% plot the normal distributions of the magnitudes derived from the scale-relationships
%%%% plot the normal distribution obatined using Mmax and its standard deviation with a black dashed line
%%%% plot Mmax with a red vertical line and its standard deviation with a dashed horizontal dotted line
clear xlabel, clear ylabel

namefig=strcat('magnitude',num2str(i)); %%name the figure to save it
figure(i)
hold on

count_pdf=1;

plot(x_range_of_mag,pdf_magnitudes(count_pdf,:),'LineWidth',1.0,'LineStyle','-','Color','b')
count_pdf=count_pdf+1;
if LAR<Length
plot(x_range_of_mag,pdf_magnitudes(count_pdf,:),'LineWidth',1.0,'LineStyle','-','Color','g')
count_pdf=count_pdf+1;
end
plot(x_range_of_mag,pdf_magnitudes(count_pdf,:),'LineWidth',1.0,'LineStyle','-','Color','r')
count_pdf=count_pdf+1;
plot(x_range_of_mag,pdf_magnitudes(count_pdf,:),'LineWidth',1.0,'LineStyle','-','Color','c')
count_pdf=count_pdf+1;
if ~isnan(mag)
plot(x_range_of_mag,pdf_magnitudes(count_pdf,:),'LineWidth',1.0,'LineStyle','-','Color','m')
end

plot(x_range_of_mag,summed_pdf_magnitudes,...
    'LineWidth',1.0,'LineStyle','--','Color',[0.3 0.3 0.3])

stem(Mmax,max(summed_pdf_magnitudes),'k','LineWidth',1.5)
line([(Mmax-sigma_Mmax) (Mmax+sigma_Mmax)], [max(summed_pdf_magnitudes) max(summed_pdf_magnitudes)],...
    'LineWidth',1.5,'LineStyle','--','Color',[0 0 0])

%% check what magnitudes entered in the graph to have the correct legend
if LAR < Length & ~isnan(mag)
  legend([ 'MMo ';'MAR ';legends_Mw;'MObs'; 'SumD'; 'Mmax'])
elseif LAR < Length & isnan(mag)
    legend([ 'MMo ';'MAR ';legends_Mw;'SumD';'Mmax'])
    elseif LAR >= Length & ~isnan(mag)
    legend([ 'MMo ';legends_Mw;'MObs'; 'SumD'; 'Mmax'])
elseif LAR >= Length & isnan(mag)
legend([ 'MMo ';legends_Mw;'SumD'; 'Mmax'])
end

xlabel('Magnitude')
ylabel('Probability density function')
hold off
title(fault)
pdf_magnitudes=[];

%%%% END OF MAIN BLOCK 2                                                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN BLOCK 3 ErrorPropagation to evaluate Tmean and alfa by errororpropagation formulas and Seismic Moment %%%%
if LAR>=Length || out_weights(2)==0; % if  LenghtFromAspectRatio >= Length given in the input or weights is zero(MAR is not used)
    L_forTmean=Length;
else
    L_forTmean=LAR;
end

Tmean=round((10.^(d+c*Mmax)./(mu*V.*L_forTmean*Width))); %average recurrence time as defined in Field,1999
Tmean =round( Tmean*(1/SCC));
Tep=[]; %average recurrence time as defined in Peruzza et al., 2010
varTep=[]; % variance of average recurrence time as defined in Peruzza et al., 2010

Tep= (10.^(d+c*Mmax)./(mu*V.*L_forTmean*Width)) +...
    ((10.^(d+c.*Mmax)*c*log(10))./(mu*V.*L_forTmean*Width)).*sigma_Mmax -...
    (10.^(d+c.*Mmax)./(mu*(V^2).*L_forTmean*Width))*dV;

Tep =Tep *(1/SCC);

varTep= (((10.^(d+c.*Mmax)*c*log(10))./(mu*V.*L_forTmean*Width))).^2 .*(sigma_Mmax.^2) +...
    ((10.^(d+c.*Mmax)./(mu*(V^2).*L_forTmean*Width))).^2*(dV^2);

varTep= varTep*(1/SCC).^2;

e=sqrt(varTep);
alfa=e./Tmean;

%%% Moment rate calculation
MomentRate(i,:)=(10.^(d+c.*((Mmax))))./Tmean;

%%%% END OF MAIN BLOCK 3                                                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN BLOCK 4 saving outputs

figname=strcat('./output_files/', outputname,'_MB_',fault);
saveas(figure(i), figname,'epsc');

out(i,:) = [i,output_Mmax_sigmaMmax(i,:),Tmean, alfa,Telap,MomentRate(i,:)];
fprintf(fidout,'%d %3.1f %3.1f %d %4.2f %d %0.4e',out(i,:));
fprintf(fidout,'%1s',blanks(1));
fprintf(fidout,'%s\n',faultname(i,:));
% temp = fprintf('%d %3.1f %3.1f %d %4.2f %d %0.4e',out(i,:))
% fprintf('%1s',blanks(1));
% fprintf('%s\n',faultname(i,:));
% MB_output = [];
% MB_output.name = faultname;
% MB_output.data = out(i,:);
% MB_output.data1 = out(:,[(1:5),12]);
% MB_output.data2 = out(:,6);
% temp = fprintf('%d %3.1f %3.1f %d %4.2f %d %0.4e',out(i,:))
% fprintf('%1s',blanks(1));
% fprintf('%s\n',faultname(i,:));
MB_output= sprintf('%d %3.1f %3.1f %d %4.2f %d %0.4e %1s %s', out(i,:),blanks(1), faultname(i,:));

    end
%%% create a LogFileMB.txt for keep trace of your calculations
LogFileMB(inputdata,outputfile,data,ScaleRel,faultname,shearmodulusfrominput,straindropfrominput,M_dM_Lengths_forLOGoutput,out,n_sigma_trunc_magnitudes,yfc)
fclose('all');
delete('./output_files/temporary_pdf_storage.txt');
fprintf('END OF CALCULATIONS\n')




