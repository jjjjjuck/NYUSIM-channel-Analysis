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

function powerSpectrum = getLosAligned(envType,powerSpectrum)
if strcmp(envType,'LOS') == true
    % Calculate the correct azimuth AoA for LOS component, which should
    % differ from its azimuth AoD by 180 degrees
    % powerSpectrum(1,4) denotes the azimuth AoD for LOS component
    clear correctAzAOA;
    if powerSpectrum(1,4) - 180 > 0
        correctAzAOA = powerSpectrum(1,4) - 180;
    else
        correctAzAOA = powerSpectrum(1,4) + 180;
    end
    % Calculate the difference between the generated azimuth AoA and
    % the correct azimuth AoA
    % powerSpectrum(1,6) is the generated azimuth AoA for LOS component
    clear dAzAOA;
    dAzAOA = powerSpectrum(1,6) - correctAzAOA;
    % Do a global shift of azimuth AoAs
    powerSpectrum(:,6) = powerSpectrum(:,6) - dAzAOA;
    powerSpectrum(:,6) = mod(powerSpectrum(:,6),360);
%     clear azAOA_temp;
%     azAOA_temp = powerSpectrum(:,6);
%     azAOA_temp(azAOA_temp < 0) = azAOA_temp(azAOA_temp < 0) + 360;
%     powerSpectrum(:,6) = azAOA_temp; 
    % Calculate the correct elevation AoA for LOS component, which
    % should be the additive inverse of the corresponding elevation AoD
    clear correctElAOA;
    correctElAOA = -powerSpectrum(1,5); 
    % Calculate the difference between the generated elevation AoA and
    % the correct elevation AoA
    % powerSpectrum(1,7) is the generated elevation AoA for LOS component
    clear dElAOA;
    dElAOA = powerSpectrum(1,7) - correctElAOA;
    % Do a global shift of elevation AoAs
    powerSpectrum(:,7) = powerSpectrum(:,7) - dElAOA;
    for i = 1:size(powerSpectrum,1)
       if  powerSpectrum(i,7) > 90
           powerSpectrum(i,7) = 180 - powerSpectrum(i,7);
       elseif powerSpectrum(i,7) < -90
           powerSpectrum(i,7) = -180 - powerSpectrum(i,7);
       end
       
    end
end