%%% NYUSIM - User License %%%

% Copyright (c) 2016-2021 New York University and NYU WIRELESS

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
function out = calPar(v1, v2, f)
% Calculate frequency-dependent channel parameters based on the measured
% parameter values at 28 GHz and 140 GHz. A linear interpolation is used 
% for frequencies between 28 GHz and 140 GHz while paramater values for 
% frequencies below 28 GHz and above 140 GHz are equal to those values at 
% 28 GHz and 140 GHz, respectively.
%
% Inputs:
%        f: carrier frequency
%        v1: parameter value at 28 GHz
%        v2: parameter value at 140 GHz
% 
% Outputs:
%        out: parameter value at input frequency f
%

getPar = @(x,y,f) f*(y-x)/(140-28)+(5*x-y)/4;
if f < 28 
    out = v1;
elseif f > 140 
    out = v2;
else
    out = getPar(v1,v2,f);
end