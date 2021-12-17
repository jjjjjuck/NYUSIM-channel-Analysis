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

function [DirRMSDelaySpread, PL_dir, PLE_dir, Pr_dir] = getDirStat(powerSpectrum,...
    theta_3dB_TX,phi_3dB_TX,theta_3dB_RX,phi_3dB_RX,TXPower,FSPL,TRDistance,d0)

%%% Obtain directional RMS delay spread, Path loss and PLE

theta_TX_d = powerSpectrum(:,4);
phi_TX_d = powerSpectrum(:,5);
theta_RX_d = powerSpectrum(:,6);
phi_RX_d = powerSpectrum(:,7);
% Number of multiapth components
nPath = size(powerSpectrum,1);
% Compute directive gains for each multipath component
PL_dir = zeros(nPath,1); 

% Directional PLE
PLE_dir = zeros(nPath,1);

% Directional RMS delay spread
DirRMSDelaySpread = zeros(nPath,1);
Pr_dir = zeros(nPath,1);
for q = 1:nPath
    % See the parameter definitions above
    [TX_Dir_Gain_Mat, RX_Dir_Gain_Mat, G_TX, G_RX] = ...
    getDirectiveGains(theta_3dB_TX,phi_3dB_TX,theta_3dB_RX,phi_3dB_RX,...
    theta_TX_d(q),phi_TX_d(q),theta_RX_d(q),phi_RX_d(q),powerSpectrum);
    [timeArray_Dir, multipathArray_Dir] = getDirPDP(powerSpectrum,...
        TX_Dir_Gain_Mat,RX_Dir_Gain_Mat);
    Pr_Lin_Dir = sum(multipathArray_Dir);
    Pr_dir(q) = 10*log10(Pr_Lin_Dir);
    meanTau = sum(timeArray_Dir.*multipathArray_Dir)/sum(multipathArray_Dir);
    meanTau_Sq=sum(timeArray_Dir.^2.*multipathArray_Dir)/sum(multipathArray_Dir);
    DirRMSDelaySpread(q) = sqrt(meanTau_Sq-meanTau^2);
    % Obtain directional path loss
    PL_dir(q) = TXPower-10*log10(Pr_Lin_Dir)+10*log10(G_TX)+10*log10(G_RX);
    % Obtain directional PLE
    PLE_dir(q) = (PL_dir(q)-FSPL)/(10*log10(TRDistance/d0));
end

end