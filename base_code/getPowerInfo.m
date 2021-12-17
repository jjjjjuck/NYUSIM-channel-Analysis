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

function [PL_dB, Pr_dBm, FSPL, PLE] = getPowerInfo(sceType,envType,f,n,SF,...
    TXPower,TRDistance,d0,p,c,u,t,RR,Pol,Fol,h_BS,folAtt,dFol)
% Apply SF and calculate received power 
if strcmp(sceType,'RMa') == false
% Note that 'getRXPower' is modified to sc_main.m
[~, PL_dB]= getRXPower(f,n,SF,TXPower,TRDistance,d0);
% RMa LOS
elseif strcmp(sceType,'RMa') == true && strcmp(envType,'LOS') == true 
    PL_dB = 20*log10(4*pi*d0*f*1e9/c) + ...
        23.1*(1-0.03*((h_BS-35)/35))*log10(TRDistance) + SF;
% RMa NLOS
elseif strcmp(sceType,'RMa') == true && strcmp(envType,'NLOS') == true 
    PL_dB = 20*log10(4*pi*d0*f*1e9/c) + ...
        30.7*(1-0.049*((h_BS-35)/35))*log10(TRDistance) + SF;
end
% Atmospheric attenuation factor
attenFactor = mpm93_forNYU(f,p,u,t,RR);
% Path loss incorporating atmospheric attenuation
PL_dB = getAtmosphericAttenuatedPL(PL_dB,attenFactor,TRDistance);
% Incorporating cross-polarization
if strcmp(Pol,'X-Pol') == true
    PL_dB = PL_dB+25;
end
% Incorporating foliage loss
if strcmp(Fol,'Yes') == true
    PL_dB = getFoliageAttenuatedPL(PL_dB,folAtt,dFol);
end      
% Calculate received power based on transmit power and path loss
Pr_dBm = TXPower - PL_dB;
% Free space path loss
FSPL = 20*log10(4*pi*d0*f*1e9/c);
% PLE
PLE = (PL_dB-FSPL)/(10*log10(TRDistance/d0));

end