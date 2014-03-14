function [metaDataH, tracesH] = ktmGetSortedOffTraces(metaData,traces)

% function [metaDataH, tracesH] = ktmGetSortedOffTraces(metaData,traces)
%
% A function to sort traces according to geophone/source abs offset 
%
% INPUT: metaData and traces
% OUTPUT:cell array offset ordered traces with {1,1}=metaData, {1,2}=traces, AND if
% opt==1, a file of traces in floats as well
%
% EXAMPLE: [metaDataH, tracesH] = ktmGetSortedOffTraces(metaData,traces)
%
% SEE ALSO: ktmGetInputMeta.m, ktmVisualInputMeta.m, ktmListInputMeta.m
%
% WJB 02/11

offset =sqrt(((double(metaData(:,24))-double(metaData(:,22)))/2).^2 + ((double(metaData(:,25))-double(metaData(:,23)))/2).^2);

[sortedH, index] = sort(offset);

metaDataH = metaData(ind,:); 
tracesH = traces(ind,:)';


if opt==1

fid=fopen(['floatTraces_hOffset',num2str(h),'.bin'],'w');

fwrite(fid,tracesH,'float');

end

