function [output] = ktmConvertXYtoOffMid(input)

% function [output] = ktmConvertXYtoOffMid(input)
%
% A function to convert btwn coord systems
%
% INPUT: [srcx, srcy, recx, recy]
%
% OUTPUT: [midx, midy, offx, offy]
%

output = [(input(:,3)+input(:,1))/2, (input(:,4)+input(:,2))/2, (input(:,3)-input(:,1))/2, (input(:,4)-input(:,2))/2];
