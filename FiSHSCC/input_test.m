% Open the file

clear; clc; close all
fid2 = fopen('./input_files/input_MB_paganica_scc.json');
% Read the file
raw = fread(fid2, inf);
% Close the file
fclose(fid2);
% Convert the data to characters
str = char(raw');
% Decode the JSON data
data = jsondecode(str);
% Access the data (example)
fault_name = fieldnames(data);
nfault = length(fault_name);
ScaleRel = [data.(fault_name{1}).ScR];
Length = [data.(fault_name{1}).Length];
Dip = [data.(fault_name{1}).Dip];
Seismogenic_Thickness = [data.(fault_name{1}).Seismogenic_Thickness];
SRmin = [data.(fault_name{1}).SRmin];
SRmax = [data.(fault_name{1}).SRmax];
Mobs = [data.(fault_name{1}).Mobs];
sdMobs = [data.(fault_name{1}).sdMobs];
Last_eq_time = [data.(fault_name{1}).Last_eq_time];
SCC = [data.(fault_name{1}).SCC];
