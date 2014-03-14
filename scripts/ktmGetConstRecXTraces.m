function [metaDataH, tracesH] = ktmGetConstRecXTraces(h,tol,opt,metaData,traces)

% function [metaDataH, tracesH] = ktmGetConstantRecXTraces(h,tol,opt,metaData,traces)
%
% A function to select traces according geophone x
%
% INPUT: recx, tol, opt, metaData and traces
% OUTPUT:cell array offsetTraces with {1,1}=metaData, {1,2}=traces, AND if
% opt==1, a file of traces and coords in floats as well
%
% EXAMPLE: [metaDataH, tracesH] = ktmGetConstantOffTraces(6500,5,1,metaData,traces)
%
% SEE ALSO: ktmGetInputMeta.m, ktmVisualInputMeta.m, ktmListInputMeta.m

bool = ((metaData(:,24)>(h-tol))&(metaData(:,24)<(h+tol))); 

[ind val]=find(bool);

%keyboard

metaDataH = metaData(ind,:); 
tracesH = traces(ind,:)';


if opt==1

fid=fopen(['../data/traces_recx',num2str(recx),'.bin'],'w');

fwrite(fid,tracesH,'float');

end

