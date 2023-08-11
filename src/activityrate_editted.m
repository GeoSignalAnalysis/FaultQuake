%%% CODE TO COMPUTE Activity rates for variuos magnitude-frequency distribution balancing Moment Rates

%%% Earthquake rates can be computed using the output file format of errorpropagation.m
%%% At 17 february we implemented 8 magnitude-frequency distribution: 
%%% 1)Single-value time independent
%%% 2)Single-value time Brownian Passage Time
%%% 3)Single-value probability given in the input file as an added last coloumn
%%% 4)Gaussian Characteristic time independent
%%% 5)Gaussian Characteristic time Brownian Passage Time
%%% 6)Gaussian Characteristic probability given in the input file as an added last coloum before name
%%% 8)ClassicalGR
%%% 7)TruncGR


global imfile
global fixed_seismicmomentrate
global Faultbehaviour
global outputname
global window
global binstep
global grinputoptions
global hazard_figures
global MB_output

%%% WARNING MESSAGES
warning off backtrace
warning off verbose
if isempty(window)
  warning('Window of observation: Where necessary you are using default value 50 years')
 window=50;
end
if isempty(binstep)
  warning('binstep: Where necessary you are using default value 0.1')
 binstep=0.1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputdata='MB_output';
Fault_behaviour=Faultbehaviour;
w=window;
bin=binstep;
hazard_fig_plot=hazard_figures;
if isempty(fixed_seismicmomentrate)
   fixed_smr=0;
else
fixed_smr=fixed_seismicmomentrate;
end

fprintf('RUNNING...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DEFINITION OF COEFFICIENTS
%%%% coefficients by Hanks & Kanamori, 1985
d=9.1; c=1.5;


%%% LOAD INPUT FILE

% in the case of models 1,2,4,5,7,8 the input file format is 
%id Mchar sdMchar Tmean alfa Telap Mo_rate(Nmyr-1) name

if Fault_behaviour == 1 |   Fault_behaviour == 2 | Fault_behaviour == 4 | Fault_behaviour == 5 | Fault_behaviour == 7 |  Fault_behaviour == 8
  
% fid = fopen(MB_output);
fid = convertCharsToStrings(MB_output);
% head1=fscanf(fid, '%s %s %s %s %s %s %s %s\n',8);
% format_string = '%d %3.1f %3.1f %d %4.2f %d %0.4e';
% data = sprintf(format_string, MB_output.data);
F=textscan(fid,'%f %f %f %f %f %f %f %s');
% fault_name = MB_output.name;
fault_name=F{8};
faultname=char(fault_name);
% nfault=size(MB_output.data,1);
% data2 = MB_output.data;
nfault=size(faultname(:,:),1);
data=[F{1},F{2},F{3},F{4},F{5},F{6},F{7},F{8}];
id=data(:,1);
mag=data(:,2);sdmag=data(:,3);
Tmean=data(:,4);alpha_val=data(:,5);Telapsed=data(:,6);
Morate_input=data(:,7);
Telapsed = (cell2mat(Telapsed));
Max_Telap=max(Telapsed); %useful for figures;
% fault_name = MB_output.name;
% nfault=size(MB_output.data,1);
% data2 = MB_output.data;

% fclose(fid);

% in the case of models 3,6 the input file format is 
%id Mchar sdMchar Tmean alfa Telap Mo_rate(Nmyr-1) Probability name
elseif Fault_behaviour == 3 |   Fault_behaviour == 6

fid = convertCharsToStrings(MB_output);
% head1=fscanf(fid, '%s %s %s %s %s %s %s %s %s\n',9);
F=textscan(fid,'%f %f %f %f %f %f %f %f %s');
fault_name=F{9};
faultname=char(fault_name);
nfault=size(faultname(:,:),1);
data=[F{1},F{2},F{3},F{4},F{5},F{6},F{7},F{8}];
id=data(:,1);
mag=data(:,2);sdmag=data(:,3);
Tmean=data(:,4);alpha_val=data(:,5);Telapsed=data(:,6);
Morate_input=data(:,7);
ProbabilityOfOccurrence=data(:,8);
Max_Telap=max(Telapsed); %useful for figures;
% fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in the case of models 7 ,8 user need a second input file with GR options:
% id, magnitude threshold and bin-spacing for magnitude-frequency distrbution
% file format is 
% id mt bin
%%% load the options for GR
ParametersGR=[];
if  ~isempty(grinputoptions);
    GRoptions=grinputoptions;
fid=fopen(GRoptions);
head1=fscanf(fid, '%s %s %s\n',3);
ParametersGRoptions=textscan(fid,'%f %f %f');
ParametersGR=[ParametersGRoptions{1},ParametersGRoptions{2},ParametersGRoptions{3}];

%ParametersGR_checkedID create a matrix sorting the GR options in the same
%way as the main input file
ParametersGR_checkedID=[];

id = (cell2mat(id));

for ck=1:length(id)
    ParametersGR_checkedID(ck,1:3) = [ParametersGR(ParametersGR(:,1)==id(ck),:)];
end
idgr=ParametersGR_checkedID(:,1);
mt=ParametersGR_checkedID(:,2);
b=ParametersGR_checkedID(:,3);
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% check FIXED TMEAN / SEISMIC MOMENT RATE

mag = (cell2mat(mag));
Tmean = (cell2mat(Tmean));
Morate_input = (cell2mat(Morate_input));

Morate_fromTmean=str2num(num2str((10.^(c.*mag+d))./Tmean,'%.4e'));
Tmean_fromMorate=round((10.^(c.*mag+d))./Morate_input);

if fixed_smr==0
    Morate=Morate_fromTmean;
elseif fixed_smr==1
    Morate=Morate_input;
end


for i=1:nfault
   if fixed_smr==0
       if Morate_fromTmean(i) ~= Morate_input(i)
    warning('Mo rate computed using M and Tmean for the fault # %d is %0.4e , different from Mo rate given in the input %0.4e',[i Morate_fromTmean(i), Morate_input(i)])
       end
   elseif fixed_smr==1
      if Tmean_fromMorate(i) ~= Tmean(i) 
    warning('Tmean computed using Mo rate for the fault # %d is %d , different from Tmean given in the input %d',[i Tmean_fromMorate(i) Tmean(i)])
       end   
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Computing of Poisson Probability and BPT probability
%%% BPT is used in models 2 and 5
%%% Depends on the Magnitude-frequency distributions these values are given
%%% as outputs
sdmag = (cell2mat(sdmag));

for i=1:nfault
    
 if Fault_behaviour >=4 
magnitude_range=(mag(i)-sdmag(i)):bin:(mag(i)+sdmag(i));
M=10.^((c.*magnitude_range)+d);
pdf_mag = pdf('Normal',magnitude_range,mag(i),sdmag(i));
total_moment=sum(pdf_mag.*M);
ratio=(Morate(i))/total_moment;
balanced_pdf_moment=ratio*pdf_mag;
CumRateMmin=sum(balanced_pdf_moment);
Tm=1/CumRateMmin;

 elseif Fault_behaviour <4 %single value models
M=10.^((c.*mag(i))+d);
Tm=1/(Morate(i)/M);
 end

% BPT hazard rate calculation
Telap=Telapsed(i);
if ~isnan(Telap)
   
if Telap>10*Tm
    Telap=10*Tm;
  warning(strcat('Telap for fault id:',blanks(1),num2str(id(i)),' is forced to be equal to 10*Tm to avoid computational problems'))
end

alpha=alpha_val(i); alpha = (cell2mat(alpha));

Hbpt_a1(1,:)=cdf('inversegaussian',(Telap+w),Tm,(Tm/(alpha^2)));
Hbpt_a2(1,:)=cdf('inversegaussian',(Telap),Tm,(Tm/(alpha^2)));
Hbpt(i)=(Hbpt_a1-Hbpt_a2)./(1-Hbpt_a2);
Hbpt(Hbpt>1)=1;
end

% POIS hazard rate calculation
Hpois(i)=1-(exp(-1*w*(1/Tm)));
Hpois(Hpois>1)=1;


%%% go here only if BPT/POISS figure is checked yes
if hazard_fig_plot == 1 & ~isnan(Telap)
   
%%%%% COMMON INPUT, x-axis of the HAZARD CURVES GRAPH
el=[];
el=(10*(Tm));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HCbpt=[];HCpois=[];
a1=[]; a2=[]; b1=[]; y1=[]; y2=[];
if ~isnan(Telap)
y1=1:1:(el+w);
y2=1:1:el;

a1(1,:)=cdf('inversegaussian',y1,Tm,(Tm/(alpha^2)));a1=a1((w+1):end);
a2(1,:)=cdf('inversegaussian',y2,Tm,(Tm/(alpha^2)));
HCbpt=(a1-a2)./(1-a2);
end
%%%% COMPUTE HAZARD FUNCTION CURVE for poisson
 HCpois= 1-(exp(-1*w*(1/Tm)));


x=y2;
figure(nfault+i);

plot(x, HCbpt,'-k') ;
hold on
line([0 max(x)], [HCpois HCpois],'LineStyle','--','Color','k') ;
plot(Telap,Hbpt(i),'or');
%plot(Tm,Hpois(i),'xb','MarkerSize',10);
line([Tm Tm], [0 max([HCbpt HCpois])],'LineStyle','-','Color','b');
hold off
xlim([0 max(el)])   
%xlim('auto')
xlabel('Telap');
ylab=strcat('P of occurrence in  ',num2str(w),' years');
ylabel(ylab);
fault=faultname(i,:);
title(fault)

if Fault_behaviour >=4 
legend('BPT','POIS','P|elap','T_M_>_=_M_m_a_x_l_b');
elseif Fault_behaviour <4 %single value models
legend('BPT','POIS','P|elap','T_M_=_M_c_e_n_t_r_a_l');
end
%%%%OUTPUT
figname=strcat('./output_files/', outputname,'_AR_Probability_',fault);
saveas(figure(nfault+i), figname,'epsc');
end    


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% choice of the magnitude-frequency distribution: call of the appropriate function

%% 1 single magnitude model, poisson
if Fault_behaviour==1
    
    [outputname]=singlePoiss(c,d,outputname,faultname,mag,Morate,id,nfault,w,Hpois);
MFD='singlePoiss';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2 single magnitude model, BPT
elseif Fault_behaviour==2
   
    [outputname]=singleBPT(outputname,mag,id,nfault,faultname,w,Hbpt);%this is the function
MFD='singleBPT';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3 single magnitude model, user defined probability
elseif Fault_behaviour==3
    
    [outputname]=singleProb(outputname,mag,id,nfault,faultname,w,ProbabilityOfOccurrence);
MFD='singleProb';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4 gaussian characteristic model, poisson
elseif Fault_behaviour==4
    
    [outputname]=CHGaussPoiss(c,d,outputname,faultname,mag,sdmag,Morate,id,nfault,w,Hpois,bin);
MFD='CHGaussPoiss';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 5 gaussian characteristic model,BPT
elseif Fault_behaviour==5
   [outputname]= CHGaussBPT(c,d,outputname,faultname,mag,sdmag,Tmean,Morate,id,nfault,w,Hbpt,bin);
MFD='CHGaussBPT';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 6 gaussian characteristic model, user defined probability
elseif Fault_behaviour==6
    
   [outputname]=CHGaussProb(c,d,outputname,faultname,mag,sdmag,Tmean,Morate,id,nfault,w,ProbabilityOfOccurrence,bin);
MFD='CHGaussProb';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 7 classical GR
elseif Fault_behaviour==7
  [outputname]=ClassicalGR(c,d,outputname,faultname,mag,mt,Morate,id,nfault,bin,b);
MFD=' ClassicalGR';   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 8 truncated GR
elseif Fault_behaviour==8
    [outputname]=TruncatedGR(c,d,outputname,faultname,mag,mt,Morate,id,nfault,bin,b);
    
MFD='TruncatedGR';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LogFileAR(inputdata,outputname,data,faultname,Fault_behaviour,MFD,window, binstep, grinputoptions,ParametersGR)

fclose('all');
fprintf('END OF CALCULATIONS\n')








