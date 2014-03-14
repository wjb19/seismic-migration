function [output] = ktmGetInputMeta(file,en,opt,range,trLen)

% function [output] = ktmGetInputMeta(file,en,opt,range,trLen)
%
% A function to strip metadata from trace files, SU/SEG-Y format
%
% INPUT: file(name) & endian & opt, where:
% opt==1, return EBCDIC reel header (if segy)
% opt==2, return binary reel header (if segy)
% opt==3, return range of trace headers (REQUIRES range to be specified)
%
% *additionally needs trLen iff *su input
%
% OUTPUT: metadata w/o labels in numerical vector/matrix
%
% EXAMPLE: out=ktmGetInputMeta('salt.c3na-b.segy','ieee-le',1)
%
% written/tested WJB 08/10 



fileID=fopen(file,'r',en);


suf = file(end-3:end);

if ((opt==1) && (suf=='segy'))

	% create lookup table, EBCDIC -> ASCII

	tbl=63*ones(1,256); 
	tbl(129:137)=['a','b','c','d','e','f','g','h','i'];
	tbl(145:153)=['j','k','l','m','n','o','p','q','r']; 
	tbl(161:169)=['~','s','t','u','v','w','x','y','z'];
	tbl(192:201)=['{','A','B','C','D','E','F','G','H','I'];
	tbl(208:217)=['}','J','K','L','M','N','O','P','Q','R']; 
	tbl(224:233)=['\',0,'S','T','U','V','W','X','Y','Z'];
	tbl(240:249)=['0','1','2','3','4','5','6','7','8','9'];
	tbl(96:97)=['-','/']; 
	tbl(64)=' '; 
	tbl(122)=':'; 
	tbl(123)='#';
	tbl(124)='@';
	

	%read f stream & write to output string  40 x 80 ASCII 		
	%characters using table lookup 

	ebcdicHDR = fread(fileID,3200,'uint8');
	ebcdicHDR=reshape(tbl(ebcdicHDR),80,40)';
	output= setstr(ebcdicHDR);

elseif ((opt==2||3) && (suf=='segy'))


	%read binary reel header, 400 bytes of ints, 2-4 bytes each
	fseek(fileID,3200,0);
	binaryReelHDR = fread(fileID,400,'uint8');
	
	%formatting

	output(1) = bin2dec(reshape(dec2bin(binaryReelHDR(1:4),8)',8*4,1)'); 
	output(2) = bin2dec(reshape(dec2bin(binaryReelHDR(5:8),8)',8*4,1)');
	output(3) = bin2dec(reshape(dec2bin(binaryReelHDR(9:12),8)',8*4,1)');

		for i=1:24
		output(i+3) = bin2dec(reshape(dec2bin(binaryReelHDR(13+(i-1)*2:14+(i-1)*2),8)',8*2,1)'); 
		end
	
end

if (opt==3 && (nargin > 3) && (range(1,2)>=range(1,1)))
	
	%in su, can't learn tr pts from line hdr
	
	if (suf(end-1:end) == 'su') && (range(1,2)-range(1,1))!=0
		
		output(8)=trLen;

		%output(10) is trace datatype; here we are simply trusting that it's int32
		%a header word would tell us explicitly but we are lazy
		output(10)=1;
	end	
	
	if (range(1,2)-range(1,1))==0
	
		output(8)=0;
		output(10)=0;
	end
		

	%data sample format length == 32 bits unless output(10)==3
 
	num2Bytes = (16 + 16*(output(10)!=3))/16; 

	%length of trace(s)

	TraceLength = num2Bytes*output(8); 
	%open file & read; should do memory check

	recordLength = (range(1,2)-range(1,1)+1);
	offset = (range(1,1)-1)*(120 + TraceLength);

	fseek(fileID,offset,3600); 
	traceHDR = fread(fileID,recordLength*(120+TraceLength),'uint16');

	traceHDR = reshape(traceHDR,(120+TraceLength),recordLength);
	traceHDR = traceHDR(1:120,:)';


	%create temp objects for formatting
	
	[m n]=size(traceHDR);

	if (en=='ieee-be')
	
	%7 4 byte numbers
	outputA = reshape(traceHDR(:,1:14),m,2,7);
	outputB = reshape((bitshift(outputA(:,1,:),16) + outputA(:,2,:)),m,7);

	%4 2 byte numbers
	%8 4 byte numbers
	outputA = reshape(traceHDR(:,19:34),m,2,8);
	outputC = reshape((bitshift(outputA(:,1,:),16) + outputA(:,2,:)),m,8);
	
	%2 2 byte numbers
	%4 4 byte numbers
	outputA = reshape(traceHDR(:,37:44),m,2,4);
	outputD = reshape((bitshift(outputA(:,1,:),16) + outputA(:,2,:)),m,4);

	else

	%7 4 byte numbers
	outputA = reshape(traceHDR(:,1:14),m,2,7);
	outputB = reshape((outputA(:,1,:) + bitshift(outputA(:,2,:),16)),m,7);

	%4 2 byte numbers
	%8 4 byte numbers
	outputA = reshape(traceHDR(:,19:34),m,2,8);
	outputC = reshape((outputA(:,1,:) + bitshift(outputA(:,2,:),16)),m,8);
	
	%2 2 byte numbers
	%4 4 byte numbers
	outputA = reshape(traceHDR(:,37:44),m,2,4);
	outputD = reshape((outputA(:,1,:) + bitshift(outputA(:,2,:),16)),m,4);

	end

	%must recast
	output = [int32(outputB) int16(traceHDR(:,15:18)) int32(outputC) int16(traceHDR(:,35:36)) int32(outputD) int16(traceHDR(:,45:n))]; 	
		

end

fclose(fileID); 









