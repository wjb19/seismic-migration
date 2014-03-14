function [Summary, finalTrace, srcX, srcY, recX, recY] = ktmWriteTraceBlockConstRecX(file,h,tol,range,opt,en,sze,imp)

% function [Summary] = ktmWriteTraceBlockConstRecX(file,h,tol,range,opt,en,sze,imp)
%
% A function to write from seismic input files, SEG-Y format, to *bin files
% (almost) ready for migration etc
%
% INPUT: file(name), recx +/- tol, desired trace range and file write mode
% (opt==1 overwrite, else append), sze is spatial block length (eg., 32),
% imp if exists is impulse test data
%
% OUTPUT: Summary of records written,  & *bin files:
% inputData.bin -> trace data uint32 (anticipate need to go ibm32 -> ieee)
% inputSrcX.bin  -> src X data float
% inputSrcY.bin  -> src Y data float
% inputRecX.bin  -> rec X data float
% inputRecY.bin  -> rec Y data float
%
% 
% EXAMPLE: ktmWriteTraceBlockConstRecX('salt.c3na-b.segy',200,10,[1,10000],1,'ieee-be',32,1)
% USES: ktmGetInputMeta, ktmGetInputTrace.m
%
% SEE ALSO: ktmListInputMeta.m, ktmVisualizeInput.m
%
% written/tested WJB 11/09 Octave 3.2.0 MacOSX 10.5.7
% updated 06/10,02/11

if opt==1
	md='w';
else
	md='a';
end

%check range

if nargin < 3 || range(1,2)<range(1,1)

error('Need to specify approp. range of traces eg., [1,10]');

end

%get the data from SEG-Y file


tic
disp('Grab Reel MetaData...')
constants= ktmGetInputMeta(file,en,2);
disp(['Done in ',num2str(toc),' s']);
tic
disp('Grab Trace MetaData...')
metaData = ktmGetInputMeta(file,en,3,range);
disp(['Done in ',num2str(toc),' s']);

tic
disp('Grab Trace Data...')
traces = ktmGetInputTrace(file,en,range);
disp(['Done in ',num2str(toc),' s']);

tic
disp(['Select Traces According to RecX = ',num2str(h),' +/- ',num2str(tol),'...'])
[metaDataH, tracesH]=ktmGetConstRecXTraces(h,tol,0,metaData,traces);

metaDataH = double(metaDataH);
RECORDS = columns(tracesH);
SAMPLES = rows(tracesH);

if (nargin ==8)


rw = ceil(SAMPLES*rand(1,RECORDS));
rc = ceil(RECORDS*rand(1,SAMPLES));


temp = zeros(size(tracesH));


for i=1:12
temp(rw(i),rc(i))=100;
end


tracesH=temp;

end


RECX = h;
TIME_STEP = constants(6);
TRACE_PTS = constants(8);
UNITS = constants(25);

Summary=[TRACE_PTS, TIME_STEP, RECORDS, RECX, UNITS];

disp('Summary: ');
disp(['TRACE_PTS: ',num2str(TRACE_PTS),' TIME_STEP: ',num2str(TIME_STEP),' RECORDS: ',num2str(RECORDS),' RECX: ',num2str(RECX)]);
disp(['MIN_SRCX: ',num2str(min(metaDataH(:,22))),' MAX_SRCX: ',num2str(max(metaDataH(:,22))),' MIN_SRCY: ',num2str(min(metaDataH(:,23))),' MAX_SRCY: ',num2str(max(metaDataH(:,23)))]);
disp(['MIN_RECX: ',num2str(min(metaDataH(:,24))),' MAX_RECX: ',num2str(max(metaDataH(:,24))),' MIN_RECY: ',num2str(min(metaDataH(:,25))),' MAX_RECY: ',num2str(max(metaDataH(:,25)))]);


if ((mod(RECORDS,sze)!=0)|| (mod(SAMPLES,sze)!=0))

newC = sze * ceil(RECORDS / sze);
newR = sze * ceil(SAMPLES / sze);

temp  = zeros(newR,newC);
tempM = zeros(newC,101);

temp(1:SAMPLES,1:RECORDS) = tracesH;
tempM(1:RECORDS,:)=metaDataH;

RECORDS=newC;
SAMPLES=newR;

traces=temp;
metaDataH=tempM;

end

finalTrace = temp;



%try and write some values

fid 	= fopen('../data/inputData.bin',md);
fid2	= fopen('../data/inputSrcX.bin',md);
fid3	= fopen('../data/inputSrcY.bin',md);
fid4	= fopen('../data/inputRecX.bin',md);
fid5	= fopen('../data/inputRecY.bin',md);
fid6    = fopen('../data/inp.txt','r');


srcX = metaDataH(:,22);
srcY = metaDataH(:,23);
recX = metaDataH(:,24);
recY = metaDataH(:,25);

if (nargin < 8)
fwrite(fid,finalTrace,'uint32');
else

fwrite(fid,finalTrace,'float');

end

fwrite(fid2,srcX,'float');
fwrite(fid3,srcY,'float');
fwrite(fid4,recX,'float');
fwrite(fid5,recY,'float');

[a b]=fscanf(fid6,"%f ",Inf);
fclose(fid6);

fid6    = fopen('../data/inp.txt','w');

%update config

a(1)=SAMPLES;
a(2)=TIME_STEP/1e6;

if opt != 1
a(3)=a(3)+RECORDS;
else
a(3)=RECORDS;
end

fprintf(fid6,"%f ",a');

fclose(fid);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
fclose(fid6);
