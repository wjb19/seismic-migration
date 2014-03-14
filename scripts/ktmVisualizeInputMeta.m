function [coord] = ktmVisualizeInputMeta(file,en,opt,range)

% function [coord] = ktmVisualizeInputMeta(file,en,opt,range)
%
% A function to visualize KTM input files, SEG-Y format
%
% INPUT: file(name) & opt, where:
% opt==1, return plot of src X/Y positions (REQUIRES range to be specified)
% opt==2, return plot of rec X/Y positions (REQUIRES range to be specified)
% opt==3, return plots of both  (REQUIRES range to be specified)
%
% OUTPUT: figure
%
% EXAMPLE: out=ktmVisualizeInputMeta('salt.c3na-b.segy','ieee-le',1,[1,10])
% USES: ktmGetInputMeta.m, ktmGetInputTrace.m
%
% SEE ALSO: ktmListInputMeta.m, ktmGetInputMeta.m, ktmGetInputTrace.m
%
% written/tested WJB 08/09 Octave 3.2.0 MacOSX 10.5.7


%check range

if nargin < 3 || range(1,2)<range(1,1)

error('Need to specify approp. range of traces eg., [1,10]');

end

if opt==1||2||3 & range(1,2)>range(1,1)

%ktmGetInputMeta(file,en,opt,range,trLen)


	temp=ktmGetInputMeta(file,en,2);
	coord = ktmGetInputMeta(file,en,3,range);

	if temp(25)==1

		Lblx=('x distance (m)');
		Lbly=('y distance (m)');

		Lbltr='(m)';

	else

		Lblx=('x distance (ft)');
		Lbly=('y distance (ft)');

		Lbltr='(ft)';

	end

end


if opt==1 & range(1,2)>range(1,1)

	plot(coord(:,22),coord(:,23),'x');
	xlabel(Lblx); ylabel(Lbly);	
	title('Source coordinates');


elseif opt==2 & range(1,2)>range(1,1)

	plot(coord(:,24),coord(:,25),'x');
	xlabel(Lblx); ylabel(Lbly);
	title('Receiver coordinates');


elseif opt==3 & range(1,2)>range(1,1)

	plot(coord(:,22),coord(:,23),'o',coord(:,24),coord(:,25),'x');
        xlabel(Lblx); ylabel(Lbly);
        title('Source (o) Receiver (x) coordinates');
	

end


