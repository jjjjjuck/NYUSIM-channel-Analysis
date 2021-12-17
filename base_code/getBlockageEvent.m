%%% NYUSIM - User License %%%

% Copyright (c) 2016-2019 New York University and NYU WIRELESS

% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the “Software”),
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software. Users shall cite 
% NYU WIRELESS publications regarding this work.

% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUTWARRANTY OF ANY KIND, EXPRESS OR 
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
% OTHER LIABILITY, WHETHER INANACTION OF CONTRACT TORT OR OTHERWISE, 
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
% OTHER DEALINGS IN THE SOFTWARE.

function [numberOfBlockage, blocksnap] = getBlockageEvent(blockChain)

% Obtain the number of blockage events, the length of each blockage
% event, the length of each state in each blockage event. 

numberOfBlockage = sum(diff(blockChain) == -3);
blocksnap = struct;
ind1 = find(diff(blockChain) == 1);
ind1 = ind1(1:3:end);
ind2 = find(diff(blockChain) == -3)+1;

for i = 1:numberOfBlockage
    indDecay = find(blockChain(ind1(i):ind2(i))==2)+ind1(i)-1;
    indShad = find(blockChain(ind1(i):ind2(i))==3)+ind1(i)-1;
    indRise = find(blockChain(ind1(i):ind2(i))==4)+ind1(i)-1;
    blocksnap.(['b',num2str(i)]).vector(1,:) = blockChain(ind1(i):ind2(i));
    blocksnap.(['b',num2str(i)]).vector(2,:) = ind1(i):ind2(i);
    blocksnap.(['b',num2str(i)]).decay = indDecay;
    blocksnap.(['b',num2str(i)]).shad = indShad;
    blocksnap.(['b',num2str(i)]).rise = indRise;
end

end

