function [freq, cumFreq, arrBound, indices_trans] = bin_var(data_in, bound)
% Find and bin the data and return each part of the data and the indices
% for these parts
% Input Args:
%   Arg 0: data_in, the data to be binned
%   Arg 1: bound, the bin width
% returns:
%   freq: the frequency of the occurance in each bin
%   cumFreq: the cumulative frequency of each occurance in a particular bin
%   arrBound: the array containing each bound
%   indices_trans: the indices of the bound data
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

