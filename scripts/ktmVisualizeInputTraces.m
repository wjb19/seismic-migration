function [data] = ktmVisualizeInputTraces(time)

% function [data] = ktmVisualizeInputTraces(time)
%
% A function to visualize trace data from ../data/inputData.bin
%
% INPUT: block & time
%
% OUTPUT: figure, data in matrix (time along columns)
%
% written/tested WJB 02/11 Octave 3.2.0 CentOS

fid = fopen('../data/inputData.bin','r');

[a b]=fread(fid,Inf,'float');

data = reshape(a,time,b/time);
