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

function [nTC,nA,nD] = getNumClusters_AOA_AOD_Indoor(mu_AOA,mu_AOD,lambda_C)
% Generate the number of time clusters, the number of AOA and AOD spatial
% lobes for the indoor scenario
%
% Inputs:
%        mu_AOA: the maximum number of AOA spatial lobes
%        mu_AOD: the maximum number of AOA spatial lobes
%        lambda_C: the mean number of time clusters
%
% Outputs:
%         nTC: the number of time clusters
%         nA: the number of AOA spatial lobes
%         nD: the number of AOD spatial lobes
%

nTC = poissrnd(lambda_C) + 1;
nA = randi(mu_AOA);
nD = randi(mu_AOD);

