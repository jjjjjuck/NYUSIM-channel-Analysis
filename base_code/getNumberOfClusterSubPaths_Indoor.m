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

function nSP = getNumberOfClusterSubPaths_Indoor(nTC,beta_s,mu_s,SPlimit)
% Generate the number of cluster subpaths for the indoor scenario
%
% Inputs:
%        nTC: the number of time clusters
%        beta_s: the beta parameter in the composite distribution - the
%                weight of the delta function
%        mu_s: the mu parameter in the composite distribution - the mean
%              value of the exponential distribution
%        SPlimit: the maximum number of subpaths in each cluster
%
% Outputs: 
%         nSP: an array of the number of subpaths per cluster
%

nSP = zeros(nTC,1);
for k = 1:nTC
    if nTC == 1
        i = 1;
    else
        i = binornd(1,beta_s);
    end
    if i == 1
        nSP(k) = floor(exprnd(mu_s)) + 1;
        while nSP(k) > SPlimit
            nSP(k) = floor(exprnd(mu_s)) + 1;
        end
    else
        nSP(k) = 1;
    end
end
% If it generates a single subpath with a single cluster, then we let it
% run again to have at least two subpaths in the generated channel.
while nTC == 1 && nSP == 1
    nSP = floor(exprnd(mu_s)) + 1;
end
