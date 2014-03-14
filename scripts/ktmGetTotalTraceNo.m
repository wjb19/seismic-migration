function [output] = ktmGetTotalTraceNo(file,en)

% function [output] = ktmGetTotalTraceNo(file,en)
%
% A function to determine total no. of traces in SEG-Y file
%
% INPUT: SEG-Y file
%
% OUTPUT: no. of traces
%
% EXAMPLE: out=ktmGetTotalTraceNo('salt.c3na-b.segy','ieee-be')
% USES: ktmGetInputMeta.m
%
% SEE ALSO: ktmListInputMeta.m, ktmVisualInputMeta.m, ktmGetInputTrace.m, ktmGetInputMeta.m
%
% written/tested WJB 08/09 Octave 3.2.0 MacOSX 10.5.7

%get trace parameters

TrParams = ktmGetInputMeta(file,en,2);

traceLength= 	(16 + 16*(TrParams(10)!=3))/8 * TrParams(8); 

%find total file length

fileID=fopen(file,'r',en);
fseek(fileID,0,'eof');

%determine no. of traces


output = (ftell(fileID)-3600) / (240+traceLength); 


fclose(fileID);


