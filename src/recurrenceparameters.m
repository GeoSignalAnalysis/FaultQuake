%%%% RECURRENCE PARAMETERS vers. 1.01
%%%% a code to compute Tmean and CV and scale-shape pairs of BPT WBL distrinutions using Monte Carlo simulation on paleoseismic dataset
%%% written by F. Visini and B. Pace with the help of L. Peruzza
%%% please refer to Peruzza et al. 2010 (JOSE) and Pace, Visini, Peruzza (2015) FiSH: Matalab tools for...
%%% when using this code or partial or modified version of it


global distinputfile
global getpath
global nsimul
global outputname
global seed_for_rnd
%%% WARNING MESSAGES
warning off backtrace
warning off verbose
if isempty(nsimul)
  warning('Number of simulations: You are using default value 10000')
 nsimul=10000;
end
%%% if exists, initialize random seed
if ~isempty(seed_for_rnd)
    seedrnd=seed_for_rnd;
   rand('seed',seedrnd);
   %rng(seedrnd);
   
else
    seedrnd=NaN;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load the paleoseismological dataset
%%% the file must be contain in each row the MIN and MAX value of the
%%% intervals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open a file for writing the output
outputfilename=strcat('output_files/',outputname, '_RP','.txt');
fidout = fopen(outputfilename,'w');


fprintf('RUNNING ...\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataset_of_sites=distinputfile;
path_of_data=getpath;
%%%load the paleoseismological dataset
%%% the files must contain in each row the MIN and MAX value of the
%%% intervals


fiddist = fopen(dataset_of_sites);
F=textscan(fiddist, '%s');
fault_name=F{1};
faultname=char(fault_name);
n_row=size(faultname(:,:),1);

%%%%%%%  LOAD input data for a site at time
for n_fault=1:n_row
    clear data, clear events, clear INT,
    clear h1,clear h2,clear h3,clear h4,
    clear k1,clear k2,clear k3,
    clear X1,clear Y1,clear X2,clear Y2,
    clear X3,clear Y3,
    clear X4
    clear param,clear parambpt,clear paramwbl,clear parampois,
    

fault=faultname(n_fault,:);
name=strcat(path_of_data,fault,'.paleo')

fid = fopen(name);
head1=fscanf(fid, '%s %s\n',2);
F=textscan(fid,'%f %f');
data=[F{1},F{2}];




%%%%%%% start with the simulation
number_of_simulations=nsimul

%%%%% create a dataset (per each simulation) of events
%for i=1:number_of_simulations
i=1;
while i < number_of_simulations
%%%%% create a dataset of events
for j=1:length(data)
events(j,:)=randi(data(j,:),1);
end 
events=sort(events,'descend');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN BLOCK 1   compute parameters                         %%%%
                                    
%%%%compute aritmetic Tmean and CV
INT=abs(diff(events));INT((INT==0))=[];
Tmean(i)=mean(INT);CV(i)=std(INT)/Tmean(i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% use a BPT distribution to compute Tmean and CV
bpt=mle(INT,'distribution','inversegaussian');
bptTmean(i)=bpt(1);
bptCV(i)=sqrt(bpt(1)/bpt(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% use a WBL distribution to compute scale and shape parameters
wbl=mle(INT,'distribution','wbl');
wblscale(i)=wbl(1);
wblshape(i)=wbl(2);

m_wbl(i)=wblscale(i).*gamma(1+(1./wblshape(i))); %mean of a weibull distribution
sd_wbl(i)=sqrt(wblscale(i).^2.*(gamma(1+2./wblshape(i))-(gamma(1+1./wblshape(i))).^2));%sdandard deviation of a weibull distribution
cv_wbl(i)=sd_wbl(i)./m_wbl(i);


%%%%%%%%%% use a POIS distribution to compute mean
pois=mle(INT,'distribution','exp');
poismean(i)=pois(1);


%%check for Inf in wblshape
if wblshape(i)==Inf
 i=i-1;
 if i==0;
     i=1;
 end
end
i=i+1;
end %%% end the simulation for a fault


%%%% END OF MAIN BLOCK 1                                         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN BLOCK 2 calculate hstograms to obtain contour graphs   %%%%

param=[CV', Tmean'];
parambpt=[bptCV',bptTmean'];
paramwbl=[wblshape',wblscale'];
paramwbl2=[cv_wbl',m_wbl'];
parampois=[poismean'];

X1=min(CV):(max(CV)-min(CV))/10:max(CV);
Y1=min(Tmean):(max(Tmean)-min(Tmean))/10:max(Tmean);
X2=min(bptCV):(max(bptCV)-min(bptCV))/10:max(bptCV);
Y2=min(bptTmean):(max(bptTmean)-min(bptTmean))/10:max(bptTmean);
X3=min(wblshape):(max(wblshape)-min(wblshape))/10:(max(wblshape));
Y3=min(wblscale):(max(wblscale)-min(wblscale))/10:max(wblscale);
X4=min(poismean):(max(poismean)-min(poismean))/10:max(poismean);
X5=min(cv_wbl):(max(cv_wbl)-min(cv_wbl))/10:max(cv_wbl);
Y5=min(m_wbl):(max(m_wbl)-min(m_wbl))/10:max(m_wbl);


n1 = hist3(param,{X1,Y1} ); % Extract histogram data; 
n1 = n1';
n2 = hist3(parambpt,{X2,Y2} ); % Extract histogram data;  
n2 = n2';
n3 = hist3(paramwbl,{X3,Y3} ); % Extract histogram data; 
n3 = n3';
n4=hist(parampois,X4);
n5 = hist3(paramwbl2,{X5,Y5} ); % Extract histogram data;  
n5 = n5';

%%%% END OF MAIN BLOCK 2                                         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN BLOCK 3 PLOT and SAVE FIGURE                           %%%%

nomefig=num2str(n_fault);
figure(n_fault)

subplot(3,2,2);contour(X1,Y1,n1);
xlabel('CV'); ylabel('Tmean');title(strcat(fault,' -arithmetic'));colorbar('East')
axis auto
subplot(3,2,4);contour(X2,Y2,n2); 
xlabel('alpha'); ylabel('Tmean');title(strcat(fault,' -BPT distribution'));colorbar('East')
axis auto
subplot(3,2,5);contour(X3,Y3,n3); 
xlabel('b'); ylabel('a');title(strcat(fault,' -WBL distribution'));colorbar('East')
axis auto
subplot(3,2,6);contour(X5,Y5,n5); 
xlabel('CV'); ylabel('Tmean');title(strcat(fault,' -WBL(Tmean,CV) distribution'));colorbar('East')
axis auto
subplot(3,2,3);hist(parampois)
xlabel('Tmean'); ylabel('count');title(strcat(fault,' -POIS distribution'))
axis auto
subplot(3,2,1);
hold on
for eq=1:size(data,1)
    line([data(eq,1) data(eq,2)], [(size(data,1)-eq+1) (size(data,1)-eq+1)],...
        'LineWidth',4,'Color',[1 0.6 0])
end
xlabel('range of occurrence (year)'); ylabel('eq#');title(strcat(fault,' -data'))
hold off
axis auto


name_of_figure=strcat('./output_files/',outputname,'_RP_',fault);
saveas(figure(n_fault),name_of_figure,'epsc')



%%%% END OF MAIN BLOCK 3                                         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN BLOCK 4 save outputs                                   %%%%

values=[];
clear ll
[h1,k1]=find(n1==max(max(n1)));
values=[X1(k1);Y1(h1)];
fprintf(fidout,'%s','ARITH: '); 
for ll=1:size(h1,1)
fprintf(fidout, 'CV= %4.2f Tmean= %3.0f',values(:,ll));
fprintf(fidout,'%1s',blanks(1));
end

values=[];
clear ll
[h2,k2]=find(n2==max(max(n2)));
values=[X2(k2);Y2(h2)];
fprintf(fidout,'%s','BPT param: '); 
for ll=1:size(h2,1)
fprintf(fidout,'alfa= %4.2f Tmean= %3.0f',values(:,ll));
fprintf(fidout,'%1s',blanks(1));
end

values=[];
clear ll
[h3,k3]=find(n3==max(max(n3)));
values=[X3(k3);Y3(h3)];
fprintf(fidout,'%s','WBL param: '); 
for ll=1:size(h3,1)
fprintf(fidout,'b= %4.2f a= %3.0f',values(:,ll));
fprintf(fidout,'%1s',blanks(1));
end

values=[];
clear ll
[h5,k5]=find(n5==max(max(n5)));
values=[X5(k5);Y5(h5)];
fprintf(fidout,'%s','WBL(CV,Tmean) param: '); 
for ll=1:size(h3,1)
fprintf(fidout,'CV= %4.2f Tmean= %3.0f',values(:,ll));
fprintf(fidout,'%1s',blanks(1));
end


values=[];
clear ll
[h4]=find(n4==max(max(n4)));
values=[X4(h4)];
fprintf(fidout,'%s','POIS param: ');  
for ll=1:size(h4,1)
fprintf(fidout,'Tmean= %3.0f',values(:,ll));
fprintf(fidout,'%1s',blanks(1));
end

fprintf(fidout,'%s\n',faultname(n_fault,:));

%%% create a LogFileRP.txt for keep trace of your calculations
LogFileRP(name,outputfilename,data,seedrnd,number_of_simulations)

end
fprintf('END OF CALCULATIONS\n')
fclose('all');

