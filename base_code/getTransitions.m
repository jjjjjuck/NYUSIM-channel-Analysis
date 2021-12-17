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

function CIR_SISO_EVO = getTransitions(nTC,numberOfSegments,co_dps,CIR_SISO_Struct,CIR_SISO_EVO)

% Realize smooth transitions between two consecutive channel segments 

for i = 1:(numberOfSegments-1)
    tc1 = nTC(i);
    sp1 = CIR_SISO_Struct.(['CIR_SISO_',num2str(i)]).numSP;
    sp1Ind = cumsum(sp1);
    
    tc2 = nTC(i+1);
    sp2 = CIR_SISO_Struct.(['CIR_SISO_',num2str(i+1)]).numSP;
    sp2Ind = cumsum(sp2);
    
    olp = max(tc1,tc2);
    olpInd = (i*co_dps-floor(olp/2)+1):(i*co_dps+ceil(olp/2)-1);
    olpInd(olpInd> length(fieldnames(CIR_SISO_EVO))) = [];
    no_olp = length(olpInd);
    trans_list = getClusterIndex(tc1,tc2);
    
    for k = 1:no_olp
        j = olpInd(k);
        trs = trans_list{k};
        no1 = sum(trs == 1);
        no2 = sum(trs == 2);
        if j > co_dps*i
            tmpCIR1 = CIR_SISO_EVO.(['Snapshot',num2str(co_dps*i)]);
            tmpCIR2 = CIR_SISO_EVO.(['Snapshot',num2str(j)]);
        else
            tmpCIR1 = CIR_SISO_EVO.(['Snapshot',num2str(j)]);
            tmpCIR2 = CIR_SISO_EVO.(['Snapshot',num2str(co_dps*i+1)]);
        end
        
        tmpCIR.pathDelays = [tmpCIR1.pathDelays(1:sp1Ind(no1));...
            tmpCIR2.pathDelays((sp2Ind(end-no2)+1):sp2Ind(end))];
        
        tmpCIR.pathPowers = [tmpCIR1.pathPowers(1:sp1Ind(no1));...
            tmpCIR2.pathPowers((sp2Ind(end-no2)+1):sp2Ind(end))];
        
        tmpCIR.pathPhases = [tmpCIR1.pathPhases(1:sp1Ind(no1));...
            tmpCIR2.pathPhases((sp2Ind(end-no2)+1):sp2Ind(end))];
        
        tmpCIR.AODs = [tmpCIR1.AODs(1:sp1Ind(no1));...
            tmpCIR2.AODs((sp2Ind(end-no2)+1):sp2Ind(end))];
        
        tmpCIR.ZODs = [tmpCIR1.ZODs(1:sp1Ind(no1));...
            tmpCIR2.ZODs((sp2Ind(end-no2)+1):sp2Ind(end))];
        
        tmpCIR.AOAs = [tmpCIR1.AOAs(1:sp1Ind(no1));...
            tmpCIR2.AOAs((sp2Ind(end-no2)+1):sp2Ind(end))];
        
        tmpCIR.ZOAs = [tmpCIR1.ZOAs(1:sp1Ind(no1));...
            tmpCIR2.ZOAs((sp2Ind(end-no2)+1):sp2Ind(end))];
        
        tmpCIR.TXPower = tmpCIR1.TXPower;
        tmpCIR.frequency = tmpCIR1.frequency;
        
        CIR_SISO_EVO.(['Snapshot',num2str(j)]) = sortCIR(tmpCIR);    
    end
end