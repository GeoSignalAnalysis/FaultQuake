function[MRLD MRA dMRLD dMRA MAR ar_coeff LAR legends_Mw]=coeff2mag(ScR,coeff,Length,Width,ARtable,mu,straindrop)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTE MAGNITUDEs FROM scale- and FROM ASPECT RATIO- relationship)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%WELLS AND COPPERSMITH 1994 EQUATIONS
% note that here Length and Width have to be expressed in km
%% check for kin of faulting

if strcmpi(ScR,'WC94-N')==1
    wc=coeff(1,:);
    ar_coeff=ARtable(1,:);
 Length_km=Length/1000;
 Width_km=Width/1000;
 
MRLD= wc(1)+ wc(2)*log10(Length_km);
MRA= wc(4) +wc(5)*log10(Length_km*Width_km);
dMRLD=wc(3);
dMRA=wc(6);

legends_Mw=['MRLD'; 'MRA '];
    
    
elseif strcmpi(ScR,'WC94-R')==1
    wc=coeff(2,:);
    ar_coeff=ARtable(2,:);
 Length_km=Length/1000;
 Width_km=Width/1000;
 
MRLD= wc(1)+ wc(2)*log10(Length_km);
MRA= wc(4) +wc(5)*log10(Length_km*Width_km);
dMRLD=wc(3);
dMRA=wc(6);

legends_Mw=['MRLD'; 'MRA '];

elseif strcmpi(ScR,'WC94-S')==1
    wc=coeff(3,:);
    ar_coeff=ARtable(3,:);
 Length_km=Length/1000;
 Width_km=Width/1000;  
    
MRLD= wc(1)+ wc(2)*log10(Length_km);
MRA= wc(4) +wc(5)*log10(Length_km*Width_km);
dMRLD=wc(3);
dMRA=wc(6);

legends_Mw=['MRLD'; 'MRA '];

elseif strcmpi(ScR,'WC94-A')==1
    wc=coeff(4,:);
    ar_coeff=ARtable(4,:);
 Length_km=Length/1000;
 Width_km=Width/1000;

MRLD= wc(1)+ wc(2)*log10(Length_km);
MRA= wc(4) +wc(5)*log10(Length_km*Width_km);
dMRLD=wc(3);
dMRA=wc(6);

legends_Mw=['MRLD'; 'MRA '];
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%LEONARD 2010 EQUATIONS


elseif strcmpi(ScR,'Le10-N')==1 | strcmpi(ScR,'Le10-R')==1 |strcmpi(ScR,'Le10-D')==1
   
    leo=coeff(1,:);
    ar_coeff=ARtable(4,:);
    
MRLDmin= (2/3)* log10(10^(leo(2)+ leo(1)*log10(Length)))-6.07;
MRAmin= (2/3)* log10(10^(leo(5) +leo(4)*log10((Length)*(Width))))-6.07;
MRLDmax= (2/3)* log10(10^(leo(3)+ leo(1)*log10(Length)))-6.07;
MRAmax= (2/3)* log10(10^(leo(6) +leo(4)*log10((Length)*(Width))))-6.07;
MRLD= MRLDmin+((MRLDmax-MRLDmin)/2);
MRA= MRAmin+((MRAmax-MRAmin)/2);
dMRLD=(MRLDmax-MRLDmin)/2;
dMRA=(MRAmax-MRAmin)/2;

legends_Mw=['MRLD'; 'MRA '];

elseif strcmpi(ScR,'Le10-S')==1
    leo=coeff(2,:);
    ar_coeff=ARtable(3,:);

MRLDmin= (2/3)* log10(10^(leo(2)+ leo(1)*log10(Length)))-6.07;
MRAmin= (2/3)* log10(10^(leo(5) +leo(4)*log10((Length)*(Width))))-6.07;
MRLDmax= (2/3)* log10(10^(leo(3)+ leo(1)*log10(Length)))-6.07;
MRAmax= (2/3)* log10(10^(leo(6) +leo(4)*log10((Length)*(Width))))-6.07;
MRLD= MRLDmin+((MRLDmax-MRLDmin)/2);
MRA= MRAmin+((MRAmax-MRAmin)/2);
dMRLD=(MRLDmax-MRLDmin)/2;
dMRA=(MRAmax-MRAmin)/2;

legends_Mw=['MRLD'; 'MRA '];

elseif strcmpi(ScR,'Le10-SCR')==1 | strcmpi(ScR,'Le10-STABLE')==1
    leo=coeff(3,:);
    ar_coeff=ARtable(4,:);
    
MRLDmin= (2/3)* log10(10^(leo(2)+ leo(1)*log10(Length)))-6.07;
MRAmin= (2/3)* log10(10^(leo(5) +leo(4)*log10((Length)*(Width))))-6.07;
MRLDmax= (2/3)* log10(10^(leo(3)+ leo(1)*log10(Length)))-6.07;
MRAmax= (2/3)* log10(10^(leo(6) +leo(4)*log10((Length)*(Width))))-6.07;
MRLD= MRLDmin+((MRLDmax-MRLDmin)/2);
MRA= MRAmin+((MRAmax-MRAmin)/2);
dMRLD=(MRLDmax-MRLDmin)/2;
dMRA=(MRAmax-MRAmin)/2;

legends_Mw=['MRLD';  'MRA '];
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%AZZARO et al. 2014 and VILLAMOR 2001 EQUATIONS (VOLCANIC CONTEXT NORMAL KIN)
% note that here Length and Width have to be expressed in km

elseif strcmpi(ScR,'Volc')==1 | strncmpi(ScR,'Az15',4)==1
    vol=coeff(1,:);
    ar_coeff=ARtable(1,:);
 Length_km=Length/1000;
 Width_km=Width/1000;
    
Mlmin= vol(1)+ vol(3)*log10(Length_km);
Mlmax= vol(2)+ vol(4)*log10(Length_km);
Mwmin= (vol(5)+ vol(6)*log10(Length_km*Width_km))-0.195;
Mwmax= (vol(5)+ vol(6)*log10(Length_km*Width_km))+0.195;
MRLD= Mlmin+((Mlmax-Mlmin)/2);%%Ml
MRA= Mwmin+((Mwmax-Mwmin)/2);%%Mw
dMRLD=(Mlmax-Mlmin)/2;%%dMl
dMRA=(Mwmax-Mwmin)/2;%%dMw

legends_Mw=['MlDA'; 'MwVi'];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aspect Ratio Control Formula (Pace and Peruzza, 2002)
% note that here Length and Width are expressed in km
Width_km=Width/1000;
LAR= ar_coeff(1)+ar_coeff(2)*(Width_km);
% check if Length from Aspect Ratio is not greater than input length

% now dimensions are in metres
%LAR=min(Length,(LAR*1000));
LAR=LAR*1000;
% calculate the moment magnitude from LAR 
MAR= (2/3)* (log10(straindrop * mu* LAR^2 *Width) -9.05);













    
    