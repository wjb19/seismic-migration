function [grid] = ktmOpenBinaryCoords()

% function [grid] = ktmOpenBinaryCoords()
%
% A function to grab and display coords from flat binary
%

fid = fopen('../data/inputRecX.bin','r');
fid1 = fopen('../data/inputRecY.bin','r');
fid2 = fopen('../data/inputSrcX.bin','r');
fid3 = fopen('../data/inputSrcY.bin','r');

[recx b]=fread(fid,Inf,'float');
[recy b]=fread(fid1,Inf,'float');
[srcx b]=fread(fid2,Inf,'float');
[srcy b]=fread(fid3,Inf,'float');

grid = [srcx, srcy, recx, recy];


