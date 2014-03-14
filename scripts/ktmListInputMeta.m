function [output] = ktmListInputMeta(file,opt,range)

% function [output] = ktmListInputMeta(file,opt,range)
%
% A function to grab metadata from KTM input files, SEG-Y format
%
% INPUT: file(name) & opt, where:
% opt==1, list EBCDIC reel header
% opt==2, list binary reel header
% opt==3, list range of trace headers (REQUIRES range to be specified)
%
% OUTPUT: metadata w/ labels in struct
%
% EXAMPLE: out=ktmListInputMeta('salt.c3na-b.segy',1)
% USES: ktmGetInputMeta.m
%
% SEE ALSO: ktmGetInputMeta.m, ktmVisualizeInput.m, ktmGetInputTrace.m, ktmIBM32toIEEE.m
%
% written/tested WJB 08/09 Octave 3.2.0 MacOSX 10.5.7




%taperCode = ['linear'; 'cosine squared'; 'other'];

%ampRecovery = ['one'; 'spherical divergence'; 'AGC'; 'other'];

%vibPolar = ['337.5 to 22.5 degrees'; '22.5 to 67.5 degrees'; '67.5 to 112.5 degrees'; '112.5 to 157.5 degrees'; '157.5 to 202.5 3degrees'; '202.5 to 247.5 degrees'; '247.5 to 292.5 degrees'; '292.5 to 337.5 degrees'];

binaryReelHDR=['Job identification number'; 'Line number'; 'Reel number'; 'Number of data traces per record'; 'Number of auxiliary traces per record'; 'Sample interval of this reel`s data in microseconds'; 'Sample interval of original field recording in microseconds'; 'Number of samples per trace for this reel`s data'; 'Number of samples per trace in original field recording'; 'Data sample format (1 = 32-bit IBM floating point, 2 = 32-bit fixed-point (integer), 3 = 16-bit fixed-point (integer), 4 = 32-bit fixed-point with gain code (integer))'; 'CDP fold (expected number of data traces per ensemble)'; 'Trace sorting code (1 = as recorded, 2 = CDP ensemble, 3 = single fold continuous profile, 4 = horizontally stacked)'; 'Vertical sum code (1 = no sum, 2  = two sum etc)'; 'Sweep frequency at start in Hertz'; 'Sweep frequency at end in Hertz'; 'Sweep length in milliseconds'; 'Sweep type code (1 = linear, 2 = parabolic, 3= exponential, 4 = other)'; 'Trace number of sweep channel'; 'Sweep trace taper length at start in milliseconds'; 'Sweep trace taper length at end in milliseconds'; 'Taper type code'; 'Correlated data traces (1 = no, 2 = yes)'; 'Binary gain recovered (1 yes 2 no)'; 'Amplitude recovery method code'; 'Measurement system (1 = meters, 2 = feet)'; 'Impulse signal polarity (increase in pressure or upward geophone case movement gives 1 (negative) or 2 (positive number))'; 'Vibratory polarity code (seismic lags pilot signal by)']; 



binaryTraceHDR=['Trace sequence number within line'; 'Trace sequence number within reel'; 'Original field record number'; 'Trace sequence number within original field record'; 'Energy source point number'; 'CDP ensemble number'; 'Trace sequence number within CDP ensemble'; 'Trace identification code (1 = seismic data, 2 = dead, 3 = dummy, 4 = time break, 5 = uphole, 6 = sweep, 7 = timing, 8 = water break, 9+ = optional use)'; 'Number of vertically summed traces yielding this trace'; 'Number of horizontally stacked traced yielding this trace'; 'Data use (1 = production, 2 = test)'; 'Distance from source point to receiver group'; 'Receiver group elevation'; 'Surface elevation at source'; 'Source depth below surface'; 'Datum elevation at receiver group'; 'Datum elevation at source'; 'Water depth at source'; 'Water depth at receiver group'; 'Scalar for elevations and depths (+ = multiplier, - = divisor)'; 'Scalar for coordinates (+ = multiplier, - = divisor)'; 'X source coordinate'; 'Y source coordinate'; 'X receiver group coordinate'; 'Y receiver group coordinate'; 'Coordinate units (1 = length in meters or feet, 2 = arc seconds)'; 'Weathering velocity'; 'Subweathering velocity'; 'Uphole time at source'; 'Uphole time at receiver group'; 'Source static correction'; 'Receiver group static correction'; 'Total static applied'; 'Lag time between end of header and time break in milliseconds'; 'Lag time between time break and shot in milliseconds'; 'Lag time beteen shot and recording start in milliseconds'; 'Start of mute time'; 'End of mute time'; 'Number of samples in this trace'; 'Sample interval of this trace in microseconds'; 'Field instrument gain type code (1 = fixed, 2 = binary, 3 = floating point, 4+ = optional use)'; 'Instrument gain constant'; 'Intrument early gain in decibels'; 'Correlated (1 = no, 2 = yes)'; 'Sweep frequency at start'; 'Sweep fequency at end'; 'Sweep length in milliseconds'; 'Sweep type code (1 = linear, 2 = parabolic, 3 = exponential,4 = other)'; 'Sweep taper trace length at start in milliseconds'; 'Sweep taper trace length at end in milliseconds'; 'Taper type code (1 = linear, 2 = cosine squared, 3 = other)'; 'Alias filter frequency'; 'Alias filter slope'; 'Notch filter frequency'; 'Notch filter slope'; 'Low cut frequency'; 'High cut frequency'; 'Low cut slope'; 'High cut slope'; 'Year data recorded'; 'Day of year'; 'Hour of day (24-hour clock)'; 'Minute of hour'; 'Second of minute'; 'Time basis (1 = local, 2 = GMT, 3 = other)'; 'Trace weighting factor for fixed-point format data'; 'Geophone group number of roll switch position one'; 'Geophone group number of first trace of original field record'; 'Geophone group number of last trace of original field record'; 'Gap size (total number of groups dropped)'; 'Overtravel associated with taper (1 = down/behind, 2 = up/ahead)']; 



if opt==2
	
	k=1;
	
	%grab numerical values from reel binary header
	reelBinaryHDR = ktmGetInputMeta(file,2);
	
	for i=1:length(reelBinaryHDR)
	
		%write struct for non-zero values

		if reelBinaryHDR(i)!=0 
		
		output{1,k} = reelBinaryHDR(i);
		output{2,k} = binaryReelHDR(i,:); 
		k++;
		end

	end

elseif opt==3 && nargin>2 && (range(1,2)>range(1,1))
	
	k=1;
	
	%grab numerical values from range of trace binary headers
	TraceBinaryHDR = ktmGetInputMeta(file,3,range);
	
	[m n]=size(TraceBinaryHDR);

	for i=1:n
	
		%write struct for non-zero values, assuming non-zero values are consistent down columns

		if TraceBinaryHDR(:,i)!=0 
		
		output{1,k} = TraceBinaryHDR(:,i);
		output{2,k} = binaryTraceHDR(i,:); 
		k++;
		end

	end

end

 

