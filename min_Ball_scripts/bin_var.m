function [freq, cumFreq, arrBound, indices_trans] = bin_var(data_in, bound)
% Find and bin the data and return each part of the data and the indices
% for these parts
% Input Args:
%   Arg 0: freq: the number of point in the bin
%   Arg 1: cumFreq: the data is required to be sorted
%   Arg 2: : the bounds that each point must fit between
% returns:
%   data_out: returns the data that fits inbetween 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the maximum value of the inside data
maxData = max(data_in);
% find the minimum value of the inside data
minData = min(data_in);

% Make an array of bounds from min to above max
arrBound = minData : bound : (maxData + bound);

diff = data_in' - arrBound;

absdiff = abs(diff);

logDiff = absdiff < bound;
indices_trans = logDiff';

freq = sum(logDiff, 1);

cumFreq = cumsum(freq, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% output_bin = cell(2, length(arrBound));
% 
% for i = 1 : length(arrBound)
%     
%     output_bin{1, i} = data_in( indices_trans(i, :) );
%     output_bin{2, i} = arrBound(i);
%     
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%