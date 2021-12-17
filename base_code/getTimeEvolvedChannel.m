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

function sc_CIR = getTimeEvolvedChannel(CIR_SISO_Struct,numberOfSegments,lengthOfSegments)

% Convert CIR struct with segments as fields to a big struct with snapshots
% as fields

sc_CIR = struct;
totInd = 1;
for segIdx = 1:numberOfSegments
  
   for j = 1:lengthOfSegments(segIdx)
       tempCIR.pathDelays = CIR_SISO_Struct.(['CIR_SISO_',num2str(segIdx)]).('Evolution').pathDelays(:,j);
       tempCIR.pathPowers = CIR_SISO_Struct.(['CIR_SISO_',num2str(segIdx)]).('Evolution').pathPowers(:,j);
       tempCIR.pathPhases = CIR_SISO_Struct.(['CIR_SISO_',num2str(segIdx)]).('Evolution').pathPhases(:,j);
       tempCIR.AODs = CIR_SISO_Struct.(['CIR_SISO_',num2str(segIdx)]).('Evolution').AODs(:,j);
       tempCIR.ZODs = CIR_SISO_Struct.(['CIR_SISO_',num2str(segIdx)]).('Evolution').ZODs(:,j);
       tempCIR.AOAs = CIR_SISO_Struct.(['CIR_SISO_',num2str(segIdx)]).('Evolution').AOAs(:,j);
       tempCIR.ZOAs = CIR_SISO_Struct.(['CIR_SISO_',num2str(segIdx)]).('Evolution').ZOAs(:,j);
       tempCIR.AOA_AOD_mapping = CIR_SISO_Struct.(['CIR_SISO_',num2str(segIdx)]).('Evolution').AOA_AOD_mapping;
       tempCIR.TXPower =  CIR_SISO_Struct.(['CIR_SISO_',num2str(segIdx)]).TXPower;
       tempCIR.frequency = CIR_SISO_Struct.(['CIR_SISO_',num2str(segIdx)]).frequency;
       sc_CIR.(['Snapshot',num2str(totInd)]) = tempCIR;
       totInd = totInd + 1;
   end
   
end