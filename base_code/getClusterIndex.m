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

function y = getClusterIndex(tc1,tc2)

r1 = ones(1,tc1)*1;
r2 = ones(1,tc2)*2;

ts = max(tc1,tc2)-1;
tmp1 = tc1;
tmp2 = tc2;
y = cell(ts,1);
for k = 1:ts
    if tmp1 > tmp2
        r1(end) = [];
        y{k} = r1;
        tmp1 = tmp1 - 1;
    elseif tmp1 < tmp2
        r1 = [r1,2];
        y{k} = r1;
        tmp2 = tmp2 - 1;
    elseif tmp1 == tmp2
        index = find(r1 == 1,1,'last');
        r1(index) = 2;
        y{k} = r1;
        tmp1 = tmp1 - 1;
        tmp2 = tmp2 - 1;
    end
end