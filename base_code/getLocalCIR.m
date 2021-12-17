%%% NYUSIM - User License %%%

% Copyright (c) 2016-2021 New York University and NYU WIRELESS

% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the �Software�),
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software. Users shall cite 
% NYU WIRELESS publications regarding this work.

% THE SOFTWARE IS PROVIDED �AS IS�, WITHOUTWARRANTY OF ANY KIND, EXPRESS OR 
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
% OTHER LIABILITY, WHETHER INANACTION OF CONTRACT TORT OR OTHERWISE, 
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
% OTHER DEALINGS IN THE SOFTWARE.

function [CIR,H,HPowers,HPhases,H_total] = getLocalCIR(CIR,TxArrayType,RxArrayType,Nt,Nr,Wt,Wr,dTxAnt,dRxAnt,RFBW) 
% This function generates the local area CIRs. TxArrayType is the type of 
% Tx antenna array (e.g., ULA, URA), RxArrayType is the type of Rx antenna 
% array, Nt is the number of Tx antennas, Nr is the number of Rx antennas, 
% Wt is # of Tx antennas in the azimuth dimention, Wr is # of Rx antennas 
% in the azimuth dimention, dTxAnt is the spacing between adjacent Tx 
% antennas in terms of wavelength, dRxAnt is the spacing between adjacent 
% Rx antennas in terms of wavelength

% Example input parameters (for test purpose only): 
% TxArrayType = 'ULA'; RxArrayType = 'URA'; Nt = 4; Nr = 10; dTxAnt = 1/2; dRxAnt = 1/2; Wt = 1; Wr = 2;

Phase = CIR.pathPhases;
n_BS = [1:1:Nt]'; n_MS = [1:1:Nr]'; 
aziAOD = CIR.AODs/360*2*pi; aziAOA = CIR.AOAs/360*2*pi;
eleAOD = CIR.ZODs/360*2*pi; eleAOA = CIR.ZOAs/360*2*pi;

nTap = size(aziAOD,1); % number of paths
Pr_lin = CIR.pathPowers; % received power in mW of each path
Pt_lin = 10^(CIR.TXPower/10);
a_BS_ULA = zeros(Nt,nTap);
a_MS_ULA = zeros(Nr,nTap);
a_BS_URA = zeros(Nt,nTap);
a_MS_URA = zeros(Nr,nTap);
for j = 1:length(n_BS)
    for k = 1:length(aziAOD)
        a_BS_ULA(j,k) = exp(1i.*(n_BS(j)-1).*2.*pi.*dTxAnt.*cos(aziAOD(k)));
    end
end
for j = 1:length(n_MS)
    for k = 1:length(aziAOA)
        a_MS_ULA(j,k) = exp(1i.*(n_MS(j)-1).*2.*pi.*dRxAnt.*cos(aziAOA(k)));
    end
end

for j = 1:length(n_BS)
    for k = 1:length(aziAOD)
        n_BS_temp = n_BS(j)-1;
        a_BS_URA(j,k) = exp(1i.*2.*pi.*dTxAnt.*...
            (mod(n_BS_temp,Wt).*cos(aziAOD(k)).*cos(eleAOD(k))+...
            fix(n_BS_temp/Wt).*sin(eleAOD(k))));
    end
end
for j = 1:length(n_MS)
    for k = 1:length(aziAOA)
        n_MS_temp = n_MS(j)-1;
        a_MS_URA(j,k) = exp(1i.*2.*pi.*dRxAnt.*...
            (mod(n_MS_temp,Wr).*cos(aziAOA(k)).*cos(eleAOA(k))+...
            fix(n_MS_temp/Wr).*sin(eleAOA(k))));
    end
end

% Calculate local area CIRs and the corresponding parameters 
H = cell(nTap,1); 
HPhases = cell(nTap,1); 
HPowers = cell(nTap,1); 
H_total = zeros(Nr,Nt);
HSmallScale = zeros(nTap,Nr);
for a = 1:nTap
    % Determine the Tx and Rx antenna array types
    if strcmp(TxArrayType,'ULA') == true && strcmp(RxArrayType,'ULA') == true 
        H{a,1} = sqrt(Pr_lin(a)/Pt_lin)*exp(1i*Phase(a))*a_MS_ULA(:,a)*a_BS_ULA(:,a)';
        H_total = H_total + H{a,1};
        % Determine the Tx and Rx antenna array types
    elseif strcmp(TxArrayType,'ULA') == true && strcmp(RxArrayType,'URA') == true 
        H{a,1} = sqrt(Pr_lin(a)/Pt_lin)*exp(1i*Phase(a))*a_MS_URA(:,a)*a_BS_ULA(:,a)';
        H_total = H_total + H{a,1};
        % Determine the Tx and Rx antenna array types
    elseif strcmp(TxArrayType,'URA') == true && strcmp(RxArrayType,'ULA') == true
        H{a,1} = sqrt(Pr_lin(a)/Pt_lin)*exp(1i*Phase(a))*a_MS_ULA(:,a)*a_BS_URA(:,a)';
        H_total = H_total + H{a,1};
        % Determine the Tx and Rx antenna array types
    elseif strcmp(TxArrayType,'URA') == true && strcmp(RxArrayType,'URA') == true
        H{a,1} = sqrt(Pr_lin(a)/Pt_lin)*exp(1i*Phase(a))*a_MS_URA(:,a)*a_BS_URA(:,a)';
        H_total = H_total + H{a,1};
    end
    % Create a SIMO channel impulse response using the first TX antenna element
    % called HSmallScale, which will be used in the Fig. 5 small-scale PDP plot
    HSmallScale(a,:) = H{a,1}(:,1); %  SIMO 
    
    CIR.HDelays{a,1} = CIR.pathDelays(a); % time delay for the ath path
    CIR.HPowers{a,1} = abs(H{a,1}).^2; % received power between each Tx antenna and each Rx antenna for the ath path
    CIR.HPhases{a,1} = CIR.pathPhases(a)+angle(H{a,1}); % phase between each Tx antenna and each Rx antenna for the ath path
    CIR.HAODs{a,1} = CIR.AODs(a); % AOD for the ath path
    CIR.HZODs{a,1} = CIR.ZODs(a); % ZOD for the ath path
    CIR.HAOAs{a,1} = CIR.AOAs(a); % AOA for the ath path
    CIR.HZOAs{a,1} = CIR.ZOAs(a); % ZOA for the ath path
    CIR.H_ensemble = H_total;
    CIR.HSmallScale = HSmallScale;
end

CIR.NumOfTxElements = Nt; % number of antenna elements in the Tx array
CIR.NumOfRxElements = Nr; % number of antenna elements in the Rx array
CIR.H = H; % H matrix
