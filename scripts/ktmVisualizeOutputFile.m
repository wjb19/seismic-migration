function [t,x,y,foo,foo2,foo3] = ktmVisualizeOutputFile(sl)

% function [output] = ktmVisualizeOutputFile(sl)
%
% A function to visualize ../data/outputData.bin using ../data/inp.txt
% assumes a file of floats, creates image planes at sl indices (x,y), (x,t), (y,t)
%
% EXAMPLE:  ktmVisualizeOutputFile([2,3,4])
%
% SEE ALSO: ktmGetInputMeta.m, ktmVisualInputMeta.m, ktmListInputMeta.m
%
% written/tested WJB 01/10,06/10 Octave 3.2.0 MacOSX 10.5.7



fid=fopen("../data/outputData.bin",'r');
[a b]=fread(fid,Inf,'float');
fclose(fid);



fid=fopen("../data/inp.txt",'r');
[c d]=fscanf(fid,"%f ",Inf);
fclose(fid);

% time/depth, y, x

output = reshape(a,c(12),c(9),c(6));
cint = 64;
colormap(gray(cint));

limt=c(12)*c(11)+c(10);
if (sl(1) > limt)
	error(['specified time slice ',num2str(sl(1)),' outside range ',num2str(limt)]);
end
t = c(10):(limt-c(10))/(c(12)-1):limt;

limx = c(6)*c(5)+c(4);
if (sl(2) > limx)
	error(['specified crossline  slice ',num2str(sl(2)),' outside range ',num2str(limx)]);
end
x = c(4):(limx-c(4))/(c(6)-1):limx;

limy = c(9)*c(8)+c(7);
if (sl(3) > limy)
	error(['specified inline  slice ',num2str(sl(3)),' outside range ',num2str(limy)]);
end
y = c(7):(limy-c(7))/(c(9)-1):limy;


if (sl(1) > 0)

	disp(['horizon at t : ',num2str(sl(1))])
	foo = squeeze(output(sl(1),:,:));
	foo = foo - min(min(foo));
	foo = cint.*foo./max(max(foo));
	image(y,x,foo);
	figure;
end

if (sl(2)>0)

	disp(['crossline at x : ',num2str(sl(2))])
	foo2 = squeeze(output(:,sl(2),:));
	foo2 = foo2 - min(min(foo2));
	foo2 = cint.*foo2./max(max(foo2));
	image(x,t,foo2);
	figure;
end

if (sl(3)>0)


	disp(['inline at y : ',num2str(sl(3))])
	foo3 = squeeze(output(:,:,sl(3)));
	foo3 = foo3 - min(min(foo3));
	foo3 = cint.*foo3./max(max(foo3));
	image(y,t,foo3);
end

