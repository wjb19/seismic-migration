function [traceData] = ktmGetInputTrace(file,en,range)

% function [traceData] = ktmGetInputTrace(file,en,range)
%
% A function to grab range of traces from KTM input files, SEG-Y format
%
% INPUT: file(name) & range of traces desired eg., [1,10]
% 
% OUTPUT: numerical traces in matrix
%
% EXAMPLE: out=ktmGetInputTrace('salt.c3na-b.segy',[1,10])
% USES: ktmGetInputMeta.m
%
% SEE ALSO: ktmListInputMeta.m
%
% written/tested WJB 08/09 Octave 3.2.0 MacOSX 10.5.7


%check range

if nargin < 2 || range(1,2)<range(1,1)

error('Need to specify approp. range of traces eg., [1,10]');

end

%grab relevant values from reel binary header; datatype etc

output = ktmGetInputMeta(file,en,2);

%grab trace values (will need to correct if not IEEE etc)

fileID=fopen(file,'r',en);
if  (range(1,2)>range(1,1))
	

	%data sample format
 
	if 	output(10)==1

		dataType = 'uint32';
		offset2  = 60;
	
	elseif 	output(10)==2

		dataType = 'int32';
		offset2  = 60;

	elseif 	output(10)==3

		dataType = 'int16';
		offset2  = 120;
	else
		dataType = 'uint32';	
		offset2  = 60;

	end
	
	dataTypeBytes = (16 + 16*(output(10)!=3))/8; 

	%length of trace(s); output(8) is no. of trace points

	TraceLength = dataTypeBytes*output(8); 

	%open file & read; should do memory check

	recordLength = (range(1,2)-range(1,1)+1);
	offset = (range(1,1)-1)*(240 + TraceLength) + 3600;

	fseek(fileID,offset,0); 
	traceData = fread(fileID,recordLength*(offset2+ output(8)),dataType);

	traceData = reshape(traceData,(offset2+output(8)),recordLength);
	traceData = traceData(offset2+1:output(8)+offset2,:)';
end



fclose(fileID);


