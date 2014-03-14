function [] = ktmWriteVelModel(vel,opt)

% function [] = ktmWriteVelModel(vel,opt)
%
% A function to write vel model to file
% 
% INPUT: velocity vector vel & opt
% (opt==1 overwrite, else append)
% 
% EXAMPLE: ktmWriteVelModel(vel)
% 
% SEE ALSO: ktmListInputMeta.m, ktmVisualizeInput.m
%
% written/tested WJB 01/10 Octave 3.2.0 MacOSX 10.5.7

if opt==1
	md='w';
else
	md='a';
end

fid 	= fopen('../data/velModel.bin',md);
fwrite(fid,vel,'float');
fclose(fid);
