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

function [powerSpectrum,hbLoss] = getHumanBlockageLoss(hbIdc,default,channel,cluster_subpath_AOAlobe_mapping,...
    cluster_subpath_AODlobe_mapping,powerSpectrum,theta_3dB_RX)

if strcmp(hbIdc,'On')
mcLen = 2e4;
t_px = 1e-3;

if strcmp(channel,'omni')
    AOA_AOD_mapping = [cluster_subpath_AOAlobe_mapping,cluster_subpath_AODlobe_mapping(:,3)];
    numAOAlobes = max(AOA_AOD_mapping(:,3));
    numAODlobes = max(AOA_AOD_mapping(:,4));
    hbLoss = zeros(numAOAlobes,numAODlobes);
    % timeArray = powerSpectrumOld(:,1); multipathArray = powerSpectrumOld(:,2); 


    for i = 1:numAOAlobes
        if isempty(find(AOA_AOD_mapping(:,3)==i,1))
            break;
        end
        allMPC_AOA = find(AOA_AOD_mapping(:,3)==i);
        maxAng = max(powerSpectrum(allMPC_AOA,6));
        minAng = min(powerSpectrum(allMPC_AOA,6));
        widthSL = maxAng-minAng; % in degree
        if widthSL > 200
            widthSL = minAng+360-maxAng;
        end
        for j = 1:numAODlobes
            if isempty(find(AOA_AOD_mapping(:,4)==j,1))
                break;
            end
            % Find the blocked MPCs
            blockMPC_AOA = find(AOA_AOD_mapping(:,3)==i);
            blockMPC_AOD = find(AOA_AOD_mapping(:,4)==j);
            blockMPC = intersect(blockMPC_AOA,blockMPC_AOD);
            if isempty(blockMPC)
                disp('No overlap');
                break;
            end
            numberOfTrace = randi(5);
            loss = zeros(numberOfTrace,mcLen);
            for m = 1:numberOfTrace 
                if strcmp(default,'Yes')

                    [mc,r] = getMarkovTrace_default(widthSL,mcLen,t_px);
                    [numberOfBlockage, blocksnap] = getBlockageEvent(mc);

                    for k = 1:numberOfBlockage

                        lengthOfDecay = length(blocksnap.(['b',num2str(k)]).decay);
                        lengthOfShad = length(blocksnap.(['b',num2str(k)]).shad);
                        lengthOfRise = length(blocksnap.(['b',num2str(k)]).rise);


                        % Decay part
                        count_decay = 1;
                        for k_decay = blocksnap.(['b',num2str(k)]).decay(1):blocksnap.(['b',num2str(k)]).decay(end)
                            loss(m,k_decay) = r*count_decay/lengthOfDecay;
                            count_decay = count_decay +1;
                        end

                        % Shadow part
                        for k_shad = blocksnap.(['b',num2str(k)]).shad(1):blocksnap.(['b',num2str(k)]).shad(end)
                            loss(m,k_shad) = r;
                        end  

                        % Rise part
                        count_rise = 1;
                        for k_rise = blocksnap.(['b',num2str(k)]).rise(1):blocksnap.(['b',num2str(k)]).rise(end)
                            loss(m,k_rise) = r*(1-count_rise/lengthOfRise);
                            count_rise = count_rise + 1;
                        end  

                    end % End of blockages for each trace


                elseif strcmp(default,'No')
                    lambdaDecay = 0.21;
                    lambdaShad = 7.88;
                    lambdaRise = 7.70;
                    lambdaUnshad = 7.67;
                    SEmean = 15.8;
                    mc = getMarkovTrace(lambdaDecay,lambdaShad,lambdaRise,lambdaUnshad,mcLen,t_px);
                    r = SEmean; % Let user input positive value

                    [numberOfBlockage, blocksnap] = getBlockageEvent(mc);

                    for k = 1:numberOfBlockage

                        lengthOfDecay = length(blocksnap.(['b',num2str(k)]).decay);
                        lengthOfShad = length(blocksnap.(['b',num2str(k)]).shad);
                        lengthOfRise = length(blocksnap.(['b',num2str(k)]).rise);

                        % Decay part
                        count_decay = 1;
                        for k_decay = blocksnap.(['b',num2str(k)]).decay(1):blocksnap.(['b',num2str(k)]).decay(end)
                            loss(m,k_decay) = r*count_decay/lengthOfDecay;
                            count_decay = count_decay +1;
                        end

                        % Shadow part
                        for k_shad = blocksnap.(['b',num2str(k)]).shad(1):blocksnap.(['b',num2str(k)]).shad(end)
                            loss(m,k_shad) = r;
                        end  

                        % Rise part
                        count_rise = 1;
                        for k_rise = blocksnap.(['b',num2str(k)]).rise(1):blocksnap.(['b',num2str(k)]).rise(end)
                            loss(m,k_rise) = r*(1-count_rise/lengthOfRise);
                            count_rise = count_rise + 1;
                        end  

                    end % End of blockages for each trace

                end % End of default on or off

            end % End of 5 blockers
            sum_loss = sum(loss,1);
            hbLoss(i,j) = sum_loss(1,randi(mcLen,1));
            linearLoss = 10^(hbLoss(i,j)/10);
            AA = exist('blockMPC');
            if AA == 0
                disp('NON exist');
            end
            powerSpectrum(blockMPC,3) = powerSpectrum(blockMPC,3)/linearLoss;

        end % End of AOD lobes
    end % End of AOA lobes
    
elseif strcmp(channel,'dir')
    [~, maxIndex] = max(powerSpectrum(:,2));
    azi_max = powerSpectrum(maxIndex,6);
    blockMPC = [];
    for impc = 1:size(powerSpectrum,1)
       if abs(powerSpectrum(impc,6)-azi_max)< theta_3dB_RX/2
           blockMPC = [blockMPC;impc];
       end
    end
    numberOfTrace = randi(5);
    loss = zeros(numberOfTrace,mcLen);
    for m = 1:numberOfTrace 
        if strcmp(default,'Yes')

            [mc,r] = getMarkovTrace_default(theta_3dB_RX,mcLen,t_px);
            [numberOfBlockage, blocksnap] = getBlockageEvent(mc);

            for k = 1:numberOfBlockage

                lengthOfDecay = length(blocksnap.(['b',num2str(k)]).decay);
                lengthOfShad = length(blocksnap.(['b',num2str(k)]).shad);
                lengthOfRise = length(blocksnap.(['b',num2str(k)]).rise);

                % Decay part
                count_decay = 1;
                for k_decay = blocksnap.(['b',num2str(k)]).decay(1):blocksnap.(['b',num2str(k)]).decay(end)
                    loss(m,k_decay) = r*count_decay/lengthOfDecay;
                    count_decay = count_decay +1;
                end

                % Shadow part
                for k_shad = blocksnap.(['b',num2str(k)]).shad(1):blocksnap.(['b',num2str(k)]).shad(end)
                    loss(m,k_shad) = r;
                end  

                % Rise part
                count_rise = 1;
                for k_rise = blocksnap.(['b',num2str(k)]).rise(1):blocksnap.(['b',num2str(k)]).rise(end)
                    loss(m,k_rise) = r*(1-count_rise/lengthOfRise);
                    count_rise = count_rise + 1;
                end  

            end % End of blockages for each trace


        elseif strcmp(default,'No')
            lambdaDecay = 0.21;
            lambdaShad = 7.88;
            lambdaRise = 7.70;
            lambdaUnshad = 7.67;
            SEmean = 15.8;
            mc = getMarkovTrace(lambdaDecay,lambdaShad,lambdaRise,lambdaUnshad,mcLen,t_px);
            r = SEmean; % Let user input positive value

            [numberOfBlockage, blocksnap] = getBlockageEvent(mc);

            for k = 1:numberOfBlockage

                lengthOfDecay = length(blocksnap.(['b',num2str(k)]).decay);
                lengthOfShad = length(blocksnap.(['b',num2str(k)]).shad);
                lengthOfRise = length(blocksnap.(['b',num2str(k)]).rise);

                % Decay part
                count_decay = 1;
                for k_decay = blocksnap.(['b',num2str(k)]).decay(1):blocksnap.(['b',num2str(k)]).decay(end)
                    loss(m,k_decay) = r*count_decay/lengthOfDecay;
                    count_decay = count_decay +1;
                end

                % Shadow part
                for k_shad = blocksnap.(['b',num2str(k)]).shad(1):blocksnap.(['b',num2str(k)]).shad(end)
                    loss(m,k_shad) = r;
                end  

                % Rise part
                count_rise = 1;
                for k_rise = blocksnap.(['b',num2str(k)]).rise(1):blocksnap.(['b',num2str(k)]).rise(end)
                    loss(m,k_rise) = r*(1-count_rise/lengthOfRise);
                    count_rise = count_rise + 1;
                end  

            end % End of blockages for each trace

        end % End of default on or off

    end % End of 5 blockers
    sum_loss = sum(loss,1);
    hbLoss = sum_loss(1,randi(mcLen,1));
    if hbLoss~= 0
        ok = 1;
    end
    linearLoss = 10^(hbLoss/10);
    powerSpectrum(blockMPC,2) = powerSpectrum(blockMPC,2)/linearLoss;

    
end % end of channel "omni" or dir
end % End of human blockage