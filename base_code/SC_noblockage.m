clear;
close all;
tic;
set(0,'defaultfigurecolor',[1 1 1])

%% User-specified basic parameters
% Set initial environment parameters
% Carrier frequency
f = 6; freq = num2str(f);

% RF bandwidth (宽带)
RFBW = 800; 

% Scenario type（新加sceType = 'InH';但不支持SC版，只有drop ）
% sceType = 'UMi';
% sceType = 'UMa';
sceType = 'UMa';
% Environment type
envType = 'NLOS'; 

% Set the initial T-R separation distance. Note that in NYUSIM v1.6.1, a
% distance range is set (dmin ~ dmax)
dmin = 100; dmax = 200; 
% dmin = 50; dmax = 50; 

% Transmit power
TXPower = 10;

% BS station height (BS is always at the origin of the simulated map)
h_BS = 10; 

% Mobile terminal height
h_MS = 1.5;

% Environment setting (pressure, humidity, foliage loss)
p = 1013.25;u = 50; t = 20;RR = 150;Fol = 'No'; dFol = 0; folAtt = 0.4; 

% TX and RX antennas setting. In the current state, only SISO is generated
% in spatial consistency procedure, but it's easy to extend to MIMO.
TxArrayType = 'ULA';Nt = 4;dTxAnt = 0.5;Wt = 1;theta_3dB_TX = 10;phi_3dB_TX = 10;
RxArrayType = 'ULA';Nr = 4;dRxAnt = 0.5;Wr = 4;theta_3dB_RX = 10;phi_3dB_RX = 10;

% type of polarization (co-polarized or cross-polarized)
Pol = 'Co-Pol';
% Pol = 'X-Pol'; % Generate total received power 会用到。

% The number of user terminal (The default value for spatial consistency procedure is 1)
N = 1; % Only one mobile
runningFolder = pwd; 
addpath(runningFolder);
%% Scenario parameters 
d0 = 1; % Free space reference distance in meters
c = physconst('LightSpeed'); % Speed of light in m/s

%% Preparation
% Structure containing generated CIRs
CIR_SISO_Struct = struct; 
CIR_MIMO_Struct = struct;

SEG_SISO_struct = struct;
SEG_MIMO_struct = struct;

% CIR_SISO_EVO = struct;
CIR_MIMO_EVO = struct;
CIR_SISO_EVO_DIR = struct;

% Set plot status
plotStatus = true; 
% Set plot rotation status
plotRotate = false; 
% Determine if spatial plot is needed 
plotSpatial = true;
FigVisibility = 'on';

scIdc = 'On';
hbIdc = 'On';

%% Set user trajectory 
% Update distance (usually set to be 1, 0.1, 0.01 m) represents the
% distance interval of two consecutive channel snapshots
% e.g. the UT moves 10 m, update distance is 1 m, then 10 channel snapshots
% are generated in total. 
d_update = 1;

% User trajectory type: 'Linear' or 'Hexagon'
% linear type is fully realized; hexagon type is almost fully realized
% except the transitions between channel segments. However, both types of
% track can be used to generate spatially correlated channel snapshots
% along the user movement without considering the transitions between
% channel segments.
trackType = 'Hexagon'; % Linear or Hexagon

% relative initial position of the UT (unit of meter) 
relPos = [0;0;h_MS]; 

% BS position (always at the origin)
TX_Pos = [0;0;h_BS];

% Moving distance (the length of user trajectory) Please use an integer
movDistance = 41; % Moving distance, unit of meter

% Moving direction where positive x-axis is 0, and positive y-axis is pi/2.
direction = 45; % Direction of the track, unit of degree

% Set the moving speed unit of meter per second
velocity = 1;

% Update time is basically update distance divided by velocity
t_update = d_update/velocity;

% The number of channel snapshots corresponding to the user movement
numberOfSnapshot = movDistance/d_update;

% Only for hexagon track
side_length = 10; % side length of the hexagon track
orient = 'Clockwise';

%%%%% An example about these distances
% A user moves 30 m. The correlation distance is 10 m, and the update
% distance is 1 m. Then, there are 3 segments (corresponding to 3 
% independent channel generations), and there are 10 snapshots in each
% segment. Thus, there are 30 snapshots in total. If the velocity is 10 m/s
% , the update time is 0.1 s, and the total moving duration is 3 s.

%% Large-scale shadow fading correlated map

% Correlation distance of SF (Please use an integer for the correlation distance, like 10, 15, 20)
d_co = 10;

% The side length of the simulated area. Since the TX (or BS) is always at
% the origin, then here the maps is from -100 to 100
area = ceil((dmax+movDistance)/100*2)*100;

%%%%% A new function of generating spatially correlated value of shadow fading
[sfMap,SF] = getSfMap(area,d_co,sceType,envType);

%% Obtain spatially correlated LOS/NLOS condition

% Usually, correlation distances for UMi, UMa, RMa are 15, 50, 60 m.
d_co_los = 15;

% Granularity of the map
% d_px = 1;

%%%%% A new function of generating spatially correlated value of LOS/NLOS
losMap = getLosMap(area,d_co_los,h_BS,h_MS,sceType);

%% The number of channel snapshots for each channel segment
co_dps = d_co/d_update;
t_dps = co_dps*t_update;

%% Channel Segments
% The moving distance may not be the multiplicity of the correlation
% distance, thus the last segment may not have full length

% Note that the envType (LOS or NLOS), sceType (UMi,UMa,RMa) do not
% change in a segment

% Generate the initial T-R separation distance
dini = getTRSep(dmin,dmax);

% Obtain the number of channel segments
numberOfSegments = ceil(movDistance/d_co);

% Obtain the vector of the lengths of channel segments 
lengthOfSegments = [co_dps*ones(1,numberOfSegments-1), mod(movDistance,d_co)/d_update];
if lengthOfSegments(end) == 0
    lengthOfSegments(end) = d_co/d_update;
end

% Generate segment-wise environment parameters
% The initial T-R separation distance for each segment (the first segment 
% has been determined.)
segDist = zeros(1,numberOfSegments); segDist(1,1) = dini;

% The environment type and scenario type of each segment (the first segment
% has been determined.)
segEnvType = cell(1,numberOfSegments); segEnvType{1,1} = envType;

% The shadow fading 
segSF = zeros(1,numberOfSegments);

% The number of time clusters in each segment (independent among segments)
nTC = zeros(1,numberOfSegments); % # of time clusters in each segment

% For loop for segments
% In each loop, an anchor channel snapshot (large-scale and small-scale 
% parameters) is first generated, then, the spatial consistency procedure
% is used to update the small-scale parameters (delays, angles, powers) of
% the rest snapshots in this segment
for segIdx = 1:numberOfSegments
    %% load channel parameters
    if strcmp(sceType,'UMi') == true && strcmp(envType,'LOS') == true 
    n = 2; SF = 4.0; mu_AOD = 1.9; mu_AOA = 1.8;X_max = 0.2;mu_tau = 123; 
    minVoidInterval = 25;sigmaCluster = 1;Gamma = 25.9; sigmaSubpath = 6; 
    gamma = 16.9; mean_ZOD = -12.6;sigma_ZOD = 5.9; std_AOD_RMSLobeAzimuthSpread = 8.5;
    std_AOD_RMSLobeElevationSpread = 2.5;distributionType_AOD = 'Gaussian'; 
    mean_ZOA = 10.8; sigma_ZOA = 5.3;std_AOA_RMSLobeAzimuthSpread = 10.5;
    std_AOA_RMSLobeElevationSpread = 11.5;distributionType_AOA = 'Laplacian';   
    % UMi NLOS
    elseif strcmp(sceType,'UMi') == true && strcmp(envType,'NLOS') == true
    n = 3.2; SF = 7.0; mu_AOD = 1.5; mu_AOA = 2.1; X_max = 0.5; mu_tau = 83;
    minVoidInterval = 25; sigmaCluster = 3; Gamma = 51.0; sigmaSubpath = 6;
    gamma = 15.5; mean_ZOD = -4.9; sigma_ZOD = 4.5; std_AOD_RMSLobeAzimuthSpread = 11.0;
    std_AOD_RMSLobeElevationSpread = 3.0; distributionType_AOD = 'Gaussian'; 
    mean_ZOA = 3.6; sigma_ZOA = 4.8; std_AOA_RMSLobeAzimuthSpread = 7.5;
    std_AOA_RMSLobeElevationSpread = 6.0; distributionType_AOA = 'Laplacian';
    % UMa LOS
    elseif strcmp(sceType,'UMa') == true && strcmp(envType,'LOS') == true 
    n = 2; SF = 4.0; mu_AOD = 1.9; mu_AOA = 1.8;X_max = 0.2; mu_tau = 123; 
    minVoidInterval = 25; sigmaCluster = 1;Gamma = 25.9; sigmaSubpath = 6; 
    gamma = 16.9; mean_ZOD = -12.6;sigma_ZOD = 5.9; std_AOD_RMSLobeAzimuthSpread = 8.5;
    std_AOD_RMSLobeElevationSpread = 2.5;distributionType_AOD = 'Gaussian'; 
    mean_ZOA = 10.8; sigma_ZOA = 5.3;std_AOA_RMSLobeAzimuthSpread = 10.5;
    std_AOA_RMSLobeElevationSpread = 11.5;distributionType_AOA = 'Laplacian'; 
    % UMa NLOS
    elseif strcmp(sceType,'UMa') == true && strcmp(envType,'NLOS') == true 
    n = 2.9; SF = 7.0; mu_AOD = 1.5; mu_AOA = 2.1; X_max = 0.5; mu_tau = 83;
    minVoidInterval = 25; sigmaCluster = 3; Gamma = 51.0; sigmaSubpath = 6;
    gamma = 15.5; mean_ZOD = -4.9; sigma_ZOD = 4.5; std_AOD_RMSLobeAzimuthSpread = 11.0;
    std_AOD_RMSLobeElevationSpread = 3.0; distributionType_AOD = 'Gaussian'; 
    mean_ZOA = 3.6; sigma_ZOA = 4.8; std_AOA_RMSLobeAzimuthSpread = 7.5;
    std_AOA_RMSLobeElevationSpread = 6.0; distributionType_AOA = 'Laplacian';
    % RMa LOS
    elseif strcmp(sceType,'RMa') == true && strcmp(envType,'LOS') == true
%     n = 2;% 这是自己加的
    SF = 1.7; mu_AOD = 1; mu_AOA = 1;X_max = 0.2; mu_tau = 123; 
    minVoidInterval = 25; sigmaCluster = 1;Gamma = 25.9; sigmaSubpath = 6; 
    gamma = 16.9; mean_ZOD = -12.6;sigma_ZOD = 5.9; std_AOD_RMSLobeAzimuthSpread = 8.5;
    std_AOD_RMSLobeElevationSpread = 2.5;distributionType_AOD = 'Gaussian'; 
    mean_ZOA = 10.8; sigma_ZOA = 5.3;std_AOA_RMSLobeAzimuthSpread = 10.5;
    std_AOA_RMSLobeElevationSpread = 11.5;distributionType_AOA = 'Laplacian';
    % RMa NLOS
    elseif strcmp(sceType,'RMa') == true && strcmp(envType,'NLOS') == true
%     n = 2;% 这是自己加的
    SF = 6.7; mu_AOD = 1; mu_AOA = 1; X_max = 0.5; mu_tau = 83;
    minVoidInterval = 25; sigmaCluster = 3; Gamma = 51.0; sigmaSubpath = 6;
    gamma = 15.5; mean_ZOD = -4.9; sigma_ZOD = 4.5; std_AOD_RMSLobeAzimuthSpread = 11.0;
    std_AOD_RMSLobeElevationSpread = 3.0; distributionType_AOD = 'Gaussian'; 
    mean_ZOA = 3.6; sigma_ZOA = 4.8; std_AOA_RMSLobeAzimuthSpread = 7.5;
    std_AOA_RMSLobeElevationSpread = 6.0; distributionType_AOA = 'Laplacian';
    end

    % Generate # of TCs,SPs,SLs
    [numberOfTimeClusters,numberOfAOALobes,numberOfAODLobes] = ...
                             getNumClusters_AOA_AOD(mu_AOA,mu_AOD,sceType);
    nTC(segIdx) = numberOfTimeClusters;
    numberOfClusterSubPaths = ...
                  getNumberOfClusterSubPaths(numberOfTimeClusters,sceType);
    nSP = numberOfClusterSubPaths;
    
    % Generate delay info
    rho_mn = getIntraClusterDelays(numberOfClusterSubPaths,X_max,sceType);
    phases_mn = getSubpathPhases(rho_mn);
    tau_n = getClusterExcessTimeDelays(mu_tau,rho_mn,minVoidInterval);
    
    % Gnerate angle info 
    % Angles between [0 360]
    [subpath_AODs, cluster_subpath_AODlobe_mapping] = ...
        getSubpathAngles(numberOfAODLobes,numberOfClusterSubPaths,mean_ZOD,...
        sigma_ZOD,std_AOD_RMSLobeElevationSpread,std_AOD_RMSLobeAzimuthSpread,...
        distributionType_AOD);
    [subpath_AOAs, cluster_subpath_AOAlobe_mapping] = ...
        getSubpathAngles(numberOfAOALobes,numberOfClusterSubPaths,mean_ZOA,...
        sigma_ZOA,std_AOA_RMSLobeElevationSpread,std_AOA_RMSLobeAzimuthSpread,...
        distributionType_AOA);
    % If it is the first channel segment running, the initial location of
    % the UT needs to be found first. 
    if segIdx == 1
        initPos = getInitPos(subpath_AODs,subpath_AOAs,segDist(segIdx),segEnvType{segIdx}, h_MS);
        % side length and orientation will not be used if the track type is
        % linear. 

        if strcmp(orient, 'Clockwise') % orientation of the hexagon track, '0'-counter;'1'-clock
            orientInd = 1;
        elseif strcmp(orient,'Counter')
            orientInd = 0;
        end
        
        [track, v_dir] = getUserTrack(trackType,initPos,movDistance,d_update,direction,side_length,orientInd);
        ss = 1+((1:numberOfSegments)-1)*co_dps;
        segInitPos = track(:,ss); 
        segDist = sqrt(sum((segInitPos-TX_Pos).^2,1));
        if sum(segDist>500) == 0
            DR = 190;
        else
            DR = 220;
        end
        Th = TXPower - DR;
        
        % Obtain LOS/NLOS condition
        for oo = 2:numberOfSegments
           losCon = losMap(round(area/2+segInitPos(2,oo)),round(area/2+segInitPos(1,oo)));
           if losCon == 1
               ttemp = 'NLOS';
           elseif losCon == 0
               ttemp = 'LOS';
           end
           segEnvType{1,oo} = ttemp; 
        end
    end
    % Obtain shadow fading from the map
    segSF(segIdx) = sfMap(round(area/2+segInitPos(2,segIdx)),round(area/2+segInitPos(1,segIdx)));
    
    % Generate total received power
    [PL_dB, Pr_dBm, FSPL, PLE] = getPowerInfo(sceType,envType,f,n,segSF(segIdx),TXPower,...
                                    segDist(segIdx),d0,p,c,u,t,RR,Pol,Fol,h_BS,folAtt,dFol); 
                                
    % Generate SP powers based on the total received power
    clusterPowers = getClusterPowers(tau_n,Pr_dBm,Gamma,sigmaCluster,Th);
    subpathPowers = ...
      getSubpathPowers(rho_mn,clusterPowers,gamma,sigmaSubpath,segEnvType(segIdx),Th);
    
    % Recover absolute timing
    t_mn = getAbsolutePropTimes(segDist(segIdx),tau_n,rho_mn);
    
    % Collect all channel information into powerSpectrum 
    % Angles between [0 360]
    powerSpectrumOld = getPowerSpectrum(numberOfClusterSubPaths,t_mn,...
                     subpathPowers,phases_mn,subpath_AODs,subpath_AOAs,Th);
    [powerSpectrum,numberOfClusterSubPaths, SubpathIndex] = ...
                                getNewPowerSpectrum(powerSpectrumOld,RFBW);
    % LOS alignment                        
    powerSpectrum = getLosAligned(envType,powerSpectrum);
    powerSpectrumOld = getLosAligned(envType,powerSpectrumOld);      
    
    % Generate a struct foe each independently generated channel segment
    CIR.pathDelays = powerSpectrumOld(:,1);
    pathPower = powerSpectrumOld(:,2);
    clear indNaN; indNaN = find(pathPower<=10^(Th/10));
    pathPower(indNaN,:) = 10^(Th/10);
    CIR.pathPowers = pathPower;
    CIR.pathPhases = powerSpectrumOld(:,3);
    CIR.AODs = powerSpectrumOld(:,4);
    CIR.ZODs = powerSpectrumOld(:,5);
    CIR.AOAs = powerSpectrumOld(:,6);
    CIR.ZOAs = powerSpectrumOld(:,7);
    CIR.frequency = freq;
    CIR.TXPower = TXPower;
    CIR.OmniPower = Pr_dBm;
    CIR.OmniPL = PL_dB;
    CIR.TRSep = segDist(segIdx);
    CIR.environment = envType;
    CIR.scenario = sceType;
    CIR.HPBW_TX = [theta_3dB_TX phi_3dB_TX];
    CIR.HPBW_RX = [theta_3dB_RX phi_3dB_RX];  
    CIR.numSP = nSP;
    % SISO CIR is stored
    CIR_SISO_Struct.(['CIR_SISO_',num2str(segIdx)]) = CIR;

    [CIR_MIMO,H,HPowers,HPhases,H_ensemble] = getLocalCIR(CIR,...
        TxArrayType,RxArrayType,Nt,Nr,Wt,Wr,dTxAnt,dRxAnt,RFBW);
    % MIMO CIR is stored
    CIR_MIMO_Struct.(['CIR_MIMO_',num2str(segIdx)]) = CIR_MIMO; 
      
    %% Time evolution
    % In a segment, the channel small-scale information will be updated for
    % each snapshot (or say for each step)
    
    % Here powerSpectrumOld for 800 MHz band width is used to ensure that the
    % number of multipath components in spatial and temporal domain are identical. 
    sc_powerSpectrum = powerSpectrumOld; 
    segTrack = track(:,(segIdx-1)*co_dps+1:(segIdx-1)*co_dps+lengthOfSegments(segIdx));
    segV = v_dir(:,(segIdx-1)*co_dps+1:(segIdx-1)*co_dps+lengthOfSegments(segIdx));
    
    no_snap = lengthOfSegments(segIdx);
    v = [cos(segV);sin(segV);zeros(1,no_snap)];
    r = zeros(3,no_snap);
    RX_Pos = segInitPos(:,segIdx);
    r(:,1) = RX_Pos - TX_Pos;
    
    % Here a geometry-based approach using multiple refleciton surfaces is
    % applied here to update angular information.
    if (strcmp(segEnvType(segIdx),'LOS'))
        no_mpc = size(sc_powerSpectrum,1)-1;
        
        % Initialize LOS and NLOS info in local coordinate system 
        sc_AOA_los = zeros(1,no_snap);sc_AOA_los(1,1) = sc_powerSpectrum(1,6);
        sc_AOD_los = zeros(1,no_snap);sc_AOD_los(1,1) = sc_powerSpectrum(1,4);
        sc_ZOA_los = zeros(1,no_snap);sc_ZOA_los(1,1) = sc_powerSpectrum(1,7);
        sc_ZOD_los = zeros(1,no_snap);sc_ZOD_los(1,1) = sc_powerSpectrum(1,5);
        sc_delay_los = zeros(1,no_snap);sc_delay_los(1,1) = sc_powerSpectrum(1,1);

        sc_AOA_nlos = zeros(no_mpc,no_snap);sc_AOA_nlos(:,1) = sc_powerSpectrum(2:end,6);
        sc_AOD_nlos = zeros(no_mpc,no_snap);sc_AOD_nlos(:,1) = sc_powerSpectrum(2:end,4);
        sc_ZOA_nlos = zeros(no_mpc,no_snap);sc_ZOA_nlos(:,1) = sc_powerSpectrum(2:end,7);
        sc_ZOD_nlos = zeros(no_mpc,no_snap);sc_ZOD_nlos(:,1) = sc_powerSpectrum(2:end,5);
        sc_delay_nlos = zeros(no_mpc,no_snap);sc_delay_nlos(:,1) = sc_powerSpectrum(2:end,1);
        
        % power is updated los and nlos together
        sc_power = zeros(no_mpc+1,no_snap);sc_power(:,1) = sc_powerSpectrum(:,2);
        sc_phase = zeros(no_mpc+1,no_snap);sc_phase(:,1) = sc_powerSpectrum(:,3);
        sc_delay = zeros(no_mpc+1,no_snap);sc_delay(:,1) = sc_powerSpectrum(:,1);
        % Initialize LOS and NLOS info in GCS
        % Convert angles into GCS
        gcs_AOD_los = zeros(1,no_snap); gcs_AOD_los(1,1) = mod(pi/2 - deg2rad(sc_AOD_los(1,1)),2*pi);
        gcs_ZOD_los = zeros(1,no_snap); gcs_ZOD_los(1,1) = pi/2 - deg2rad(sc_ZOD_los(1,1));
        gcs_AOA_los = zeros(1,no_snap); gcs_AOA_los(1,1) = mod(pi/2 - deg2rad(sc_AOA_los(1,1)),2*pi);
        gcs_ZOA_los = zeros(1,no_snap); gcs_ZOA_los(1,1) = pi/2 - deg2rad(sc_ZOA_los(1,1));

        gcs_AOD_nlos = zeros(no_mpc,no_snap); gcs_AOD_nlos(:,1) = mod(pi/2 - deg2rad(sc_AOD_nlos(:,1)),2*pi);
        gcs_ZOD_nlos = zeros(no_mpc,no_snap); gcs_ZOD_nlos(:,1) = pi/2 - deg2rad(sc_ZOD_nlos(:,1));
        gcs_AOA_nlos = zeros(no_mpc,no_snap); gcs_AOA_nlos(:,1) = mod(pi/2 - deg2rad(sc_AOA_nlos(:,1)),2*pi);
        gcs_ZOA_nlos = zeros(no_mpc,no_snap); gcs_ZOA_nlos(:,1) = pi/2 - deg2rad(sc_ZOA_nlos(:,1));
 
        xBern = randi(2,no_mpc,1);
        for t = 2:no_snap
            
            % Update LOS component
            r(:,t) = r(:,t-1) + v(:,t-1)*t_update;
            sc_delay_los(1,t) = norm(r(:,t))/c*1e9;
            
            % delta is the difference term between two snapshots 
            deltaAOD_los = v(:,t-1)'*[-sin(gcs_AOD_los(1,t-1));cos(gcs_AOD_los(1,t-1));0]*t_update;
            gcs_AOD_los(1,t) = gcs_AOD_los(1,t-1) + deltaAOD_los/(c*sc_delay_los(1,t-1)*1e-9*sin(gcs_ZOD_los(1,t-1)));

            deltaZOD_los = v(:,t-1)'*[cos(gcs_ZOD_los(1,t-1))*cos(gcs_AOD_los(1,t-1));cos(gcs_ZOD_los(1,t-1))*sin(gcs_AOD_los(1,t-1));-sin(gcs_ZOD_los(1,t-1))]*t_update;
            gcs_ZOD_los(1,t) = gcs_ZOD_los(1,t-1) + deltaZOD_los/(c*sc_delay_los(1,t-1)*1e-9);

            deltaAOA_los = v(:,t-1)'*[-sin(gcs_AOA_los(1,t-1));cos(gcs_AOA_los(1,t-1));0]*t_update;
            gcs_AOA_los(1,t) = gcs_AOA_los(1,t-1) - deltaAOA_los/(c*sc_delay_los(1,t-1)*1e-9*sin(gcs_ZOA_los(1,t-1)));

            deltaZOA_los = v(:,t-1)'*[cos(gcs_ZOA_los(1,t-1))*cos(gcs_AOA_los(1,t-1));cos(gcs_ZOA_los(1,t-1))*sin(gcs_AOA_los(1,t-1));-sin(gcs_ZOA_los(1,t-1))]*t_update;
            gcs_ZOA_los(1,t) = gcs_ZOA_los(1,t-1) + deltaZOA_los/(c*sc_delay_los(1,t-1)*1e-9);
            
            % Update delay
            azi = gcs_AOA_nlos(:,t-1);
            ele = gcs_ZOA_nlos(:,t-1);
            r_hat = [cos(ele).*sin(azi), cos(ele).*cos(azi), sin(ele)];
            deltaDist = r_hat*v(:,t-1)*t_update;
            sc_delay_nlos(:,t) = sc_delay_nlos(:,t-1)-deltaDist/c*1e9;            
            [sc_tau_n,sc_rho_mn] = getDelayInfo([sc_delay_los(t);sc_delay_nlos(:,t)],nSP);
            
            % Update power
            sfc = sfMap(round(area/2+track(1,t+(segIdx-1)*co_dps)),round(area/2+track(2,t+(segIdx-1)*co_dps)));
            dc = sqrt(sum(track(:,t+(segIdx-1)*co_dps).^2));
            [~, Prc,~,~] = getPowerInfo(sceType,envType,f,n,sfc,TXPower,...
                                    dc,d0,p,c,u,t,RR,Pol,Fol,h_BS,folAtt,dFol);
            sc_power_cluster = getClusterPowers(sc_tau_n,Prc,Gamma,sigmaCluster,Th);
            power_temp = getSubpathPowers(sc_rho_mn,sc_power_cluster,gamma,sigmaSubpath,segEnvType(segIdx),Th);
            sc_power(:,t) = structToList(power_temp,nSP);
            
            % Update phase
            sc_delay(:,t) = [sc_delay_los(1,t);sc_delay_nlos(:,t)];
            dt_delay = sc_delay(:,t)-sc_delay(:,t-1);
            sc_phase(:,t) = mod(sc_phase(:,t-1) + dt_delay*2*pi*f*1e-3, 2*pi);
            
            % Update NLOS components one by one
            for i_path = 1:no_mpc
            
                tempBern = xBern(i_path);
                deltaRS = gcs_AOA_nlos(i_path,t-1)+(-1)^tempBern*gcs_AOD_nlos(i_path,t-1)+tempBern*pi;
                v_RS = mod(deltaRS+(-1)^tempBern*v(1:2,t-1),2*pi);
                v_RS = [v_RS;0];
                
                deltaAOD = v_RS'*[-sin(gcs_AOD_nlos(i_path,t-1));cos(gcs_AOD_nlos(i_path,t-1));0]*t_update;
                gcs_AOD_nlos(i_path,t) = gcs_AOD_nlos(i_path,t-1) + deltaAOD/(c*sc_delay_nlos(i_path,t-1)*1e-9*sin(gcs_ZOD_nlos(i_path,t-1)));

                deltaZOD = v_RS'*[cos(gcs_ZOD_nlos(i_path,t-1))*cos(gcs_AOD_nlos(i_path,t-1));cos(gcs_ZOD_nlos(i_path,t-1))*sin(gcs_AOD_nlos(i_path,t-1));-sin(gcs_ZOD_nlos(i_path,t-1))]*t_update;
                gcs_ZOD_nlos(i_path,t) = gcs_ZOD_nlos(i_path,t-1) + deltaZOD/(c*sc_delay_nlos(i_path,t-1)*1e-9);

                deltaAOA = v_RS'*[-sin(gcs_AOA_nlos(i_path,t-1));cos(gcs_AOA_nlos(i_path,t-1));0]*t_update;
                gcs_AOA_nlos(i_path,t) = gcs_AOA_nlos(i_path,t-1) - deltaAOA/(c*sc_delay_nlos(i_path,t-1)*1e-9*sin(gcs_ZOA_nlos(i_path,t-1)));

                deltaZOA = v_RS'*[cos(gcs_ZOA_nlos(i_path,t-1))*cos(gcs_AOA_nlos(i_path,t-1));cos(gcs_ZOA_nlos(i_path,t-1))*sin(gcs_AOA_nlos(i_path,t-1));-sin(gcs_ZOA_nlos(i_path,t-1))]*t_update;
                gcs_ZOA_nlos(i_path,t) = gcs_ZOA_nlos(i_path,t-1) + deltaZOA/(c*sc_delay_nlos(i_path,t-1)*1e-9);
                
            end
           
        end
        
        % Change angles back to local coordinate system
        sc_AOD_los = mod(rad2deg(pi/2 - gcs_AOD_los),360);sc_AOD_nlos = mod(rad2deg(pi/2 - gcs_AOD_nlos),360);
        sc_ZOD_los = rad2deg(pi/2 - gcs_ZOD_los);sc_ZOD_nlos = rad2deg(pi/2 - gcs_ZOD_nlos);
        sc_AOA_los = mod(rad2deg(pi/2 - gcs_AOA_los),360);sc_AOA_nlos = mod(rad2deg(pi/2 - gcs_AOA_nlos),360);
        sc_ZOA_los = rad2deg(pi/2 - gcs_ZOA_los);sc_ZOA_nlos = rad2deg(pi/2 - gcs_ZOA_nlos);
        
        % Save updated info of all snapshots
        evoCIR.pathDelays = [sc_delay_los;sc_delay_nlos];
        sc_pathPower = sc_power;
        % Deal with low powers
        for tt = 1:no_snap
            clear indNaN; 
            indNaN = find(sc_pathPower(:,tt)<=10^(Th/10));
            sc_pathPower(indNaN,tt) = 10^(Th/10);
        end
        evoCIR.pathPowers = sc_pathPower;
        evoCIR.pathPhases = sc_phase;
        evoCIR.AODs = [sc_AOD_los;sc_AOD_nlos];
        evoCIR.ZODs = [sc_ZOD_los;sc_ZOD_nlos];
        evoCIR.AOAs = [sc_AOA_los;sc_AOA_nlos];
        evoCIR.ZOAs = [sc_ZOA_los;sc_ZOA_nlos];
        evoCIR.no_snap = no_snap;
        evoCIR.AOA_AOD_mapping = [cluster_subpath_AODlobe_mapping cluster_subpath_AOAlobe_mapping(:,3)];
        CIR_SISO_Struct.(['CIR_SISO_',num2str(segIdx)]).('Evolution') = evoCIR;
        
        
    elseif (strcmp(segEnvType(segIdx),'NLOS'))
        no_mpc = size(sc_powerSpectrum,1);
        sc_AOA_nlos = zeros(no_mpc,no_snap);sc_AOA_nlos(:,1) = sc_powerSpectrum(:,6);
        sc_AOD_nlos = zeros(no_mpc,no_snap);sc_AOD_nlos(:,1) = sc_powerSpectrum(:,4);
        sc_ZOA_nlos = zeros(no_mpc,no_snap);sc_ZOA_nlos(:,1) = sc_powerSpectrum(:,7);
        sc_ZOD_nlos = zeros(no_mpc,no_snap);sc_ZOD_nlos(:,1) = sc_powerSpectrum(:,5);
        sc_delay_nlos = zeros(no_mpc,no_snap);sc_delay_nlos(:,1) = sc_powerSpectrum(:,1);
        sc_power = zeros(no_mpc,no_snap);sc_power(:,1) = sc_powerSpectrum(:,2);  
        sc_phase = zeros(no_mpc,no_snap);sc_phase(:,1) = sc_powerSpectrum(:,3);
        
        gcs_AOD_nlos = zeros(no_mpc,no_snap); gcs_AOD_nlos(:,1) = mod(pi/2 - deg2rad(sc_AOD_nlos(:,1)),2*pi);
        gcs_ZOD_nlos = zeros(no_mpc,no_snap); gcs_ZOD_nlos(:,1) = pi/2 - deg2rad(sc_ZOD_nlos(:,1));
        gcs_AOA_nlos = zeros(no_mpc,no_snap); gcs_AOA_nlos(:,1) = mod(pi/2 - deg2rad(sc_AOA_nlos(:,1)),2*pi);
        gcs_ZOA_nlos = zeros(no_mpc,no_snap); gcs_ZOA_nlos(:,1) = pi/2 - deg2rad(sc_ZOA_nlos(:,1));
        
        xBern = randi(2,no_mpc,1);
        
        for t = 2:no_snap

            % Update delay
            azi = gcs_AOA_nlos(:,t-1);
            ele = gcs_ZOA_nlos(:,t-1);
            r_hat = [cos(ele).*sin(azi), cos(ele).*cos(azi), sin(ele)];
            deltaDist = r_hat*v(:,t-1)*t_update;
            sc_delay_nlos(:,t) = sc_delay_nlos(:,t-1)-deltaDist/c*1e9;            
%             [sc_tau_n,sc_rho_mn] = getDelayInfo([sc_delay_los(t);sc_delay_nlos(:,t)],nSP);
            [sc_tau_n,sc_rho_mn] = getDelayInfo(sc_delay_nlos(:,t),nSP);


            % Update power
            sfc = sfMap(round(area/2+track(1,t+(segIdx-1)*co_dps)),round(area/2+track(2,t+(segIdx-1)*co_dps)));
            dc = sqrt(sum(track(:,t+(segIdx-1)*co_dps).^2));
            [~, Prc,~,~] = getPowerInfo(sceType,envType,f,n,sfc,TXPower,...
                                    dc,d0,p,c,u,t,RR,Pol,Fol,h_BS,folAtt,dFol);
            sc_power_cluster = getClusterPowers(sc_tau_n,Pr_dBm,Gamma,sigmaCluster,Th);
            power_temp = getSubpathPowers(sc_rho_mn,sc_power_cluster,gamma,sigmaSubpath,segEnvType(segIdx),Th);
            sc_power(:,t) = structToList(power_temp,nSP);
            
            % Update phase
            dt_delay = sc_delay_nlos(:,t)-sc_delay_nlos(:,t-1);
            sc_phase(:,t) = mod(sc_phase(:,t-1) + dt_delay*2*pi*f*1e-3, 2*pi);
            
            % Update angles (NLOS components)
            for i_path = 1:no_mpc

                tempBern = xBern(i_path);
                deltaRS = gcs_AOA_nlos(i_path,t-1)+(-1)^tempBern*gcs_AOD_nlos(i_path,t-1)+tempBern*pi;
                v_RS = mod(deltaRS+(-1)^tempBern*v(1:2,t-1),2*pi);
                v_RS = [v_RS;0];
                
                deltaAOD = v_RS'*[-sin(gcs_AOD_nlos(i_path,t-1));cos(gcs_AOD_nlos(i_path,t-1));0]*t_update;
                gcs_AOD_nlos(i_path,t) = gcs_AOD_nlos(i_path,t-1) + deltaAOD/(c*sc_delay_nlos(i_path,t-1)*1e-9*sin(gcs_ZOD_nlos(i_path,t-1)));

                deltaZOD = v_RS'*[cos(gcs_ZOD_nlos(i_path,t-1))*cos(gcs_AOD_nlos(i_path,t-1));cos(gcs_ZOD_nlos(i_path,t-1))*sin(gcs_AOD_nlos(i_path,t-1));-sin(gcs_ZOD_nlos(i_path,t-1))]*t_update;
                gcs_ZOD_nlos(i_path,t) = gcs_ZOD_nlos(i_path,t-1) + deltaZOD/(c*sc_delay_nlos(i_path,t-1)*1e-9);

                deltaAOA = v_RS'*[-sin(gcs_AOA_nlos(i_path,t-1));cos(gcs_AOA_nlos(i_path,t-1));0]*t_update;
                gcs_AOA_nlos(i_path,t) = gcs_AOA_nlos(i_path,t-1) - deltaAOA/(c*sc_delay_nlos(i_path,t-1)*1e-9*sin(gcs_ZOA_nlos(i_path,t-1)));

                deltaZOA = v_RS'*[cos(gcs_ZOA_nlos(i_path,t-1))*cos(gcs_AOA_nlos(i_path,t-1));cos(gcs_ZOA_nlos(i_path,t-1))*sin(gcs_AOA_nlos(i_path,t-1));-sin(gcs_ZOA_nlos(i_path,t-1))]*t_update;
                gcs_ZOA_nlos(i_path,t) = gcs_ZOA_nlos(i_path,t-1) + deltaZOA/(c*sc_delay_nlos(i_path,t-1)*1e-9);

            end 
        end
        
        % Change angles back to local coordinate system
        sc_AOD_nlos = mod(rad2deg(pi/2 - gcs_AOD_nlos),360);
        sc_ZOD_nlos = rad2deg(pi/2 - gcs_ZOD_nlos);
        sc_AOA_nlos = mod(rad2deg(pi/2 - gcs_AOA_nlos),360);
        sc_ZOA_nlos = rad2deg(pi/2 - gcs_ZOA_nlos);
        
        evoCIR.pathDelays = sc_delay_nlos;
        % Deal with low powers
        for tt = 1:no_snap
            clear indNaN; 
            indNaN = find(sc_power(:,tt)<=10^(Th/10));
            sc_power(indNaN,tt) = 10^(Th/10);
        end
        evoCIR.pathPowers = sc_power;
        evoCIR.pathPhases = sc_phase;
        evoCIR.AODs = sc_AOD_nlos;
        evoCIR.ZODs = sc_ZOD_nlos;
        evoCIR.AOAs = sc_AOA_nlos;
        evoCIR.ZOAs = sc_ZOA_nlos;
        evoCIR.no_snap = no_snap;
        evoCIR.AOA_AOD_mapping = [cluster_subpath_AOAlobe_mapping cluster_subpath_AODlobe_mapping(:,3)];
        CIR_SISO_Struct.(['CIR_SISO_',num2str(segIdx)]).('Evolution') = evoCIR;
        
    end
     
end % End of a channel segment

%% Put all snapshots together
CIR_SISO_EVO = getTimeEvolvedChannel(CIR_SISO_Struct,numberOfSegments,lengthOfSegments);

%% Smooth transitions
CIR_SISO_EVO = getTransitions(nTC,numberOfSegments,co_dps,CIR_SISO_Struct,CIR_SISO_EVO);

TR3D = sqrt(sum(track.^2));
TR2D = sqrt(sum(track(1:2,:).^2));
omniPL = zeros(numberOfSnapshot,1);
omniPr = zeros(numberOfSnapshot,1);
omniDS = zeros(numberOfSnapshot,1);
KFactor = zeros(numberOfSnapshot,1);

for i_snap = 1:numberOfSnapshot
    
    % Omnidirectional channels
    sdf = sfMap(round(area/2+track(1,i_snap)),round(area/2+track(2,i_snap)));
    [PL_dB, Pr_dBm, FSPL, PLE] = getPowerInfo(sceType,envType,f,n,segSF(segIdx),TXPower,...
                                    segDist(segIdx),d0,p,c,u,t,RR,Pol,Fol,h_BS,folAtt,dFol); 
    omniPL(i_snap,1) = PL_dB;
    omniPr(i_snap,1) = Pr_dBm;
    CIR_tmp = CIR_SISO_EVO.(['Snapshot',num2str(i_snap)]);
    multipathArray = CIR_tmp.pathPowers;
    Pr = 10*log10(multipathArray);
    xmaxInd = find(Pr>Th);
    Pr = Pr(xmaxInd);
    timeArray = CIR_tmp.pathDelays;
    timeArray = timeArray(xmaxInd);
    multipathArray = multipathArray(xmaxInd);
    meanTau = sum(timeArray.*multipathArray)/sum(multipathArray);
    meanTau_Sq = sum(timeArray.^2.*multipathArray)/sum(multipathArray);
    RMSDelaySpread = sqrt(meanTau_Sq-meanTau^2);
    omniDS(i_snap) = RMSDelaySpread;
    KFactor(i_snap) = 10*log10(max(multipathArray)/(sum(multipathArray)-max(multipathArray)));
    
end

%% MIMO CIR for each channel snapshot
for i_snap = 1:numberOfSnapshot
    CIR_tmp = CIR_SISO_EVO.(['Snapshot',num2str(i_snap)]);
    [CIR_MIMO_tmp,~,~,~,~] = getLocalCIR(CIR_tmp,...
        TxArrayType,RxArrayType,Nt,Nr,Wt,Wr,dTxAnt,dRxAnt,RFBW);
    % MIMO CIR is stored
    CIR_MIMO_EVO.(['Snapshot',num2str(i_snap)]) = CIR_MIMO_tmp;
end

%% Directional CIR for each channel snapshot
DirPDPInfo = [];
for i_snap = 1:numberOfSnapshot
    CIR_tmp = CIR_SISO_EVO.(['Snapshot',num2str(i_snap)]);
    ps = [CIR_tmp.pathDelays, CIR_tmp.pathPowers, CIR_tmp.pathPhases,...
        CIR_tmp.AODs, CIR_tmp.ZODs, CIR_tmp.AOAs, CIR_tmp.ZOAs];
    TRd = sqrt(sum(track(1:2,i_snap).^2));
    
    [DirRMSDelaySpread, PL_dir, PLE_dir, Pr_dir] = getDirStat(ps,...
    theta_3dB_TX,phi_3dB_TX,theta_3dB_RX,phi_3dB_RX,TXPower,FSPL,TRd,d0);
    
    % Plot the strongest directional PDP
    [maxP, maxIndex] = max(ps(:,2));
    
    % Angles for use
    theta_TX_d = CIR_tmp.AODs;
    phi_TX_d = CIR_tmp.ZODs;
    theta_RX_d = CIR_tmp.AOAs;
    phi_RX_d = CIR_tmp.ZOAs;
    
    % Get directive antenna gains
    [TX_Dir_Gain_Mat, RX_Dir_Gain_Mat, G_TX, G_RX] = getDirectiveGains(theta_3dB_TX,...
        phi_3dB_TX,theta_3dB_RX,phi_3dB_RX,theta_TX_d(maxIndex),phi_TX_d(maxIndex),...
        theta_RX_d(maxIndex),phi_RX_d(maxIndex),ps);

    % Recover the directional PDP
    [timeArray_Dir, multipathArray_Dir] = getDirPDP(ps,...
        TX_Dir_Gain_Mat,RX_Dir_Gain_Mat);
    
%     Pr_Dir = 10^(sum(multipathArray_Dir)/10);
    meanTau = sum(timeArray_Dir.*multipathArray_Dir)/sum(multipathArray_Dir);
    meanTau_Sq = sum(timeArray_Dir.^2.*multipathArray_Dir)/sum(multipathArray_Dir);
    rmsDS_best = sqrt(meanTau_Sq-meanTau^2);
    DirPDP = [timeArray_Dir, multipathArray_Dir];
    
    CIR_tmp.pathPowers_BestDir = multipathArray_Dir; % Best direction
    CIR_tmp.rmsDS = DirRMSDelaySpread;
    CIR_tmp.PL_dir = PL_dir;
    CIR_tmp.PLE_dir = PLE_dir;
    CIR_tmp.rmsDS_BestDir = rmsDS_best;
    CIR_tmp.Pr_dir = Pr_dir;
    CIR_SISO_EVO_DIR.(['Snapshot',num2str(i_snap)]) = CIR_tmp;
    
    onefill = ones(length(Pr_dir),1);
    DirPDPInfo_temp = [i_snap*onefill,TR3D(i_snap)*onefill,Pr_dir,PL_dir,DirRMSDelaySpread];
    DirPDPInfo = vertcat(DirPDPInfo,DirPDPInfo_temp);
end


toc;