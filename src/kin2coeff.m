function[coeff,ARtable]=kin2coeff(ScR)

if strncmpi(ScR,'WC94',4)==1
%%%%% coefficients by Wells & Coppersmith, 1994
%%%%% coefficients to calculate  Magnitude from length and area
%% matrix, row 1 for normal; row 2 for reverse; row 3 for strikeslip; row 4 for ALL
%% each row contains aRLD, bRLD, sdRLD, aRA, bRA, sdRA
coeff=[4.34, 1.54, 0.31, 3.93, 1.02, 0.25;
       4.49, 1.49, 0.26, 4.33, 0.90, 0.25;
       4.33, 1.49, 0.24, 3.98, 1.02, 0.23;
       4.38, 1.49, 0.26, 4.07, 0.98, 0.24];
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif strncmpi(ScR,'Le10',4)==1
%%%%% coefficients by Leonard, 2010
%%%%% coefficients to calculate  Moment from length and area
%% matrix, row 1 for DS dip slip; row 2 for SS strike slip; row 3 for SCR Stable Continental Regions
%% each row contains aRLD, bRLDmin, bRLDmax, aRA, bRAmin, bRAmax
%% RLD=L RA=A (see Table 5 in Leonard, 2010)
%% DS: A>5,500m2 SS:A45,000 m2 SCR: A>2,500 m2        
    coeff=[2.5, 7.53, 8.51, 1.5, 5.69, 6.6;
           1.5, 12.01, 12.88, 1.5, 5.69, 6.47;
           2.5, 7.87, 8.28, 1.5, 6.22, 6.52];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif strncmpi(ScR,'Az15',4)==1 | strncmpi(ScR,'Volc',4)==1
%%%%% coefficients by D'Amico and Azzaro, 2014 and Villamor 2001
%%%%% coefficients to calculate Ml (D'Amico) and Mw (Villamor) from
%%%%% surface rupture length (D'Amico) and RA (Villamor)
%% matrix, row 1 for D'Amico; row 2 for Villamor;
%% each row contains aSRLmin, aSRLmax, bSRLmin, bSRLmax,
%% SRL=surface rupture lenght        
    coeff=[3.239, 3.543, 1.662, 2.49, 3.39, 1.33];
    
    
end



%%%% ASPECT RATIO
%%coefficients by Pace et al.,2002 (BGTA)
%% matrix, row 1 for normal; row 2 for reverse; row 3 for strikeslip; row 4 for ALL DIP
%% each row contains aAS, bAS, sdAS
%% these coefficients are ALWAYS used for computing Magnitudes
ARtable=[3.0939,  1.2501, 0.25;
    -4.4543, 2.1992, 0.25;
    -7.096,  2.9807, 0.25;
    -2.3725, 1.9354, 0.25];