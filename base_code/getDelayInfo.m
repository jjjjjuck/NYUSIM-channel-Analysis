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

function [tau_n,rho_mn] = getDelayInfo(DelayList,nSP)
    % Note that DelayList is LOS + NLOS
    % get excess delay
    nTC = size(nSP,1);
    tau = DelayList'-min(DelayList);
    for i = 1:nTC
        tmp(1,i) = sum(nSP(1:i));
    end 
    firstComponentIdx = [0 tmp(1:(nTC-1))] + 1;
    tau_n = tau(firstComponentIdx);
    
    edgeIdx = [firstComponentIdx, size(DelayList,1)+1];
    rho_mn = struct;
    for cIdx=1:nTC
        sp = tau(edgeIdx(cIdx):(edgeIdx(cIdx+1)-1));
        spDelay = sp - min(sp);
        str = ['c',num2str(cIdx)];
        rho_mn.(str) = spDelay;
    end
end