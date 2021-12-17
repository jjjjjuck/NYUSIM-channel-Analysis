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

function [sfMap,SF] = getSfMap(area,d_co,sceType,envType)

% Obtain the map of spatially correlated shadow fading 

% The below values are abstracted from the field measurements. 
if strcmp(sceType,'UMi') == true && strcmp(envType,'LOS') == true 
    SF = 4.0;
elseif strcmp(sceType,'UMi') == true && strcmp(envType,'NLOS') == true
    SF = 7.0;
elseif strcmp(sceType,'UMa') == true && strcmp(envType,'LOS') == true 
    SF = 4.0;
elseif strcmp(sceType,'UMa') == true && strcmp(envType,'NLOS') == true 
    SF = 7.0;
elseif strcmp(sceType,'RMa') == true && strcmp(envType,'LOS') == true
    SF = 1.7;
elseif strcmp(sceType,'RMa') == true && strcmp(envType,'NLOS') == true
    SF = 6.7;
end

% Generate grid
N = area + 1;
d_px = 1;
delta = d_co/d_px;

% The filter coefficient is considered as 0 when the distance is beyond
% 4*d_co.
M = 8*delta+1;
h = zeros(M);
InitMap = randn(N+M-1);

% Generate the filter
for i = 1:M
    for j = 1:M
        h(i,j) = exp(-sqrt(((M+1)/2-i).^2+((M+1)/2-j).^2)/d_co);
    end
end

% The mean of the shadow fading
mu = 0;

% Filter the intial map
CorrMapPad = conv2(InitMap,h,'same');

% Corp the map back to the original size
CorrMap = CorrMapPad((M+1)/2:(M+1)/2+N-1,(M+1)/2:(M+1)/2+N-1);

% Do normalization
% Calculate the actual mean
mu0 = mean(CorrMap(:));
% Calculate the actual variance
sigma0 = sqrt(var(CorrMap(:)));
% Scale the Correlated Map 
sfMap = CorrMap*(SF/sigma0)+(mu-mu0);

end


