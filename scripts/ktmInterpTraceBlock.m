function [interpTrace] = ktmInterpTraceBlock(binTrace,N,opt)

% function [interpTrace] = ktmInterpTraceBlock(binTrace,N,opt)
%
% A function to interp trace blocks according to 2^N grid
%
% INPUT: binTrace, N and opt, where:
% opt==1 interp space
% opt==2 interp time
%
% binTrace is a cell array with binTrace{1,1}= [metaData1; metaData2; ..], 
% binTrace{1,2}=[trace1; trace2; ...]
% OUTPUT: same cell array with interp temporal points (opt==2) OR for opt==1
% interpTrace{1,1} = [metaData1; metaData2; ..], interpTrace{1,2}= new X grid, interpTrace{1,3}=new Y grid, 
% interpTrace{1,4} = new Traces
% 
% EXAMPLE: interpTrace = ktmInterpTraceBlock(trblock,1,2)
% 
%
% SEE ALSO: ktmListInputMeta.m, ktmVisualizeInput.m, ktmGetInputTrace.m, ktmIBM32toIEEE.m
%
% written/tested WJB 08/09 Octave 3.2.0 MacOSX 10.5.7


coord = binTrace{1,1};

traces = binTrace{1,2};


[m n]=size(coord);

[M NN]=size(traces);




if opt==1

% need to learn ordering of data to produce appropriate grid
% tol == tolerance in original spatial sampling; perfect grid == 0 

tol = 0.1; rowLength=0; colLength=0;

if coord(2,25)-coord(1,25) < tol
	
		
	%row (x) major format; learn length of rows
		
	rowLength=1; valDiff = tol; 

	while (abs(valDiff)<= tol) && (rowLength < m)
	
	valDiff = (coord(rowLength+1,25) - coord(rowLength,25));

		if abs(valDiff) <= tol

		rowLength++;

		end

	end

	%calculate no. of traces we can take to interp, must be floor(m/rowLength)* rowLength

	colLength = floor(m/rowLength);

	X = reshape(coord(1:floor(m/rowLength)* rowLength,24),rowLength,colLength);
	Y = reshape(coord(1:floor(m/rowLength)* rowLength,25),rowLength,colLength);

	tracesStack = reshape(traces(1:floor(m/rowLength)* rowLength,:),rowLength,colLength,NN);

	
	
	
else

	%column (y) major format; learn length of columns
		
	colLength=1; valDiff = tol; 

	while (abs(valDiff)<= tol) && (colLength < m)
	
	valDiff = (coord(rowLength+1,24) - coord(rowLength,24));

		if abs(valDiff) <= tol

		colLength++;

		end

	end

	%calculate no. of traces we can take to interp, must be floor(m/colLength)* colLength

	rowLength = floor(m/colLength);


	X = reshape(coord(1:floor(m/rowLength)* rowLength,24),rowLength,colLength);
	Y = reshape(coord(1:floor(m/rowLength)* rowLength,25),rowLength,colLength);

	tracesStack = reshape(traces(1:floor(m/rowLength)* rowLength,:),rowLength,colLength,NN);
end

%interp according to X/Y




	UpBndX = max(coord(:,24)); LoBndX = min(coord(:,24));
	UpBndY = max(coord(:,25)); LoBndY = min(coord(:,25));


	[XX YY]=meshgrid(LoBndX:(UpBndX-LoBndX)/(2^N-1):UpBndX, LoBndY:(UpBndY-LoBndY)/(2^N-1):UpBndY);


	%now step through spatially separated traces, interp to new mesh

	for i=1:NN

		newTraces(:,:,i) = interp2(X,Y,tracesStack(:,:,i),XX',YY');

	end


	interpTrace{1,1}= binTrace{1,1}; 
	interpTrace{1,2}=XX';
	interpTrace{1,3}=YY';
	interpTrace{1,4}=newTraces;

	

%temperal interp; for ease of recalculating dwell time, strictly work with 1/2^N

else



	origTime = 0:NN-1;

	newTime =  0:1/(2^N):NN-1;

	
	
	%now interp traces, to new intervals

	for i=1:M

		newTraces(i,:) = interp1(origTime,traces(i,:),newTime);

	end
	
	interpTrace{1,1}= binTrace{1,1}; 
	interpTrace{1,2}=newTraces;


end
