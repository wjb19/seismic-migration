function [metaDataH, tracesH] = ktmGetConstOffTraces(h,tol,opt,metaData,traces)

% function [metaDataH, tracesH] = ktmGetConstantOffTraces(h,tol,opt,metaData,traces)
%
% A function to select traces according to particular geophone/source offset h+/-tol 
%
% INPUT: h, tol, opt, metaData and traces
% OUTPUT:cell array offsetTraces with {1,1}=metaData, {1,2}=traces, AND if
% opt==1, a file of traces in floats as well
%
% EXAMPLE: [metaDataH, tracesH] = ktmGetConstantOffTraces(-140,5,1,metaData,traces)
%
% SEE ALSO: ktmGetInputMeta.m, ktmVisualInputMeta.m, ktmListInputMeta.m
%
% written/tested WJB 08/09 Octave 3.2.0 MacOSX 10.5.7
% updated 06/10 not returning cell array anymore

offset =round(sqrt(((metaData(:,24)-metaData(:,22))/2).^2 + ((metaData(:,25)-metaData(:,23))/2).^2));


bool = ((offset>(h-tol))&(offset<(h+tol))); 

[ind val]=find(bool);

%keyboard

metaDataH = metaData(ind,:); 
tracesH = traces(ind,:)';


if opt==1

fid=fopen(['floatTraces_hOffset',num2str(h),'.bin'],'w');

fwrite(fid,tracesH,'float');

end

