function [rhoTrace] = ktmApplyRhoFilter(trace,dw,opt)

% function [rhoTrace] = ktmApplyRhoFilter(trace,dw,opt)
%
% A function to apply the rho filter to trace data block
%
% INPUT: time domain traces & dw (dwell time in micro s) & opt, where:
% opt==1, apply 2D
% opt==2, apply 2.5D
% pot==3, apply 3D
% 
% OUTPUT: filtered data in time domain
%
% EXAMPLE: rhoTrace=ktmApplyRhoFilter(traces,8e3,1)
%
% SEE ALSO: ktmGetInputMeta.m, ktmVisualInputMeta.m, ktmListInputMeta.m
% REF: "From the Hagedoorn imaging technique to Kirchoff migration and inversion"
% N. Bleistein, S. H. Gray, Geophysical Prospecting 2001 (49) 629-643
%
% written/tested WJB 08/09 Octave 3.2.0 MacOSX 10.5.7

[n m]=size(trace);

%FFT is fast on 2^n, pad to M

M = 2^ceil(log(m)/log(2));


%Nyquist thm, frequencies

Dr=dw*1e-6;

omega=repmat((2*pi*(-M/2:M/2-1)./Dr),n,1)';

%FFT; put zero freq in middle of spec (fftshift), take FFT of transpose (fft goes down columns)

fTrace = fftshift(fft(trace',M)); 

if opt==1

	%do 2D filter

	rhoTrace = real(ifft(fftshift(abs(omega).*fTrace)))';
	rhoTrace = rhoTrace(:,1:m);

elseif opt==2

	%do 2.5D filter

	rhoTrace = real(ifft(fftshift(sqrt(abs(omega)).*exp((i*pi/4).*sign(omega)).*fTrace)))';
	rhoTrace = rhoTrace(:,1:m);

else

	%do 3D filter

	rhoTrace = real(ifft(fftshift(i.*omega.*fTrace)))';
	rhoTrace = rhoTrace(:,1:m);
	

	
end

 
