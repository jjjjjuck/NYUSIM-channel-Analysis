clear;
close all;
tic;%与toc结合使用记录程序运行时间
set(0,'defaultfigurecolor',[1 1 1])

%% User-specified basic parameters 用户基础参数
% Set initial environment parameters
% Carrier frequency 频率
f = 28; freq = num2str(f); % 单位GHz

% RF bandwidth 带宽  
RFBW = 500; % 单位MHz 

% Scenario type 环境
sceType = 'UMi';

% Environment type 环境类型
envType = 'NLOS'; 

% Set the initial T-R separation distance 收发距离
dmin = 100; dmax = 200;  % 单位m
% dmin = 50; dmax = 50; 

% Transmit power 发送功率
TXPower = 10; % 单位dBm

% BS station height 基站高度 (BS is always at the origin of the simulated map)
h_BS = 10; 

% Mobile terminal height 终端高度
h_MS = 1.5;

% Environment setting (pressure, humidity, foliage loss) 环境设置 
p = 1013.25;u = 50; t = 20;RR = 150;Fol = 'No'; dFol = 0; folAtt = 0.4; % do these used?

% TX and RX antennas setting. 天线设置
TxArrayType = 'ULA';Nt = 16;dTxAnt = 0.5;Wt = 1;theta_3dB_TX = 10;phi_3dB_TX = 10;
RxArrayType = 'ULA';Nr = 16;dRxAnt = 0.5;Wr = 1;theta_3dB_RX = 10;phi_3dB_RX = 10;

% TxArrayType = 'ULA';Nt = 128;dTxAnt = 0.5;Wt = 1;theta_3dB_TX = 10;phi_3dB_TX = 10;
% RxArrayType = 'ULA';Nr = 64;dRxAnt = 0.5;Wr = 1;theta_3dB_RX = 10;phi_3dB_RX = 10;

% type of polarization (co-polarized or cross-polarized)损耗设置
Pol = 'Co-Pol';
% Pol = 'X-Pol';

% The number of user terminal 用户终端数量 (The default value for spatial consistency procedure is 1)
N = 1; % Only one mobile 目前只支持单链路

%% Scenario parameters 环境参数
d0 = 1; % Free space reference distance in meters 自由空间参考距离 
c = 3e8; % Speed of light in m/s 光速

%% Preparation 预定义变量
% Structure containing generated CIRs 存储CIR
CIR_SISO_Struct = struct; 
CIR_MIMO_Struct = struct;

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

%% Set user trajectory 用户移动设定
% Update distance 两次通道快照间的距离，为1或0.1m
d_update = 0.1;% 最小好像是0.0001

% 移动单位为每秒/m
velocity = 1;

% Moving distance 用户移动的距离，为整数
movDistance = 20; % Moving distance, unit of meter

%更新时间
t_update = d_update/velocity;

% User trajectory type: 'Linear'线性型 or 'Hexagon'六角形
trackType = 'Linear'; % Linear or Hexagon

% relative initial position of the UT， 用户终端的相对初始位置(unit of meter) 
relPos = [0;0;h_MS]; 

% BS position (always at the origin)基站的相对初始位置
TX_Pos = [0;0;h_BS];

% Moving direction where positive x-axis is 0, and positive y-axis is pi/2.
% x轴正方向为0度，y轴正方向为90度
direction = 45; % Direction of the track, unit of degree

% Only for hexagon track
side_length = 10; % side length of the hexagon track
orient = 'Clockwise';

% The number of channel snapshots corresponding to the user movement
%通道快照的总数
numberOfSnapshot = movDistance/d_update;

%%%%% An example about these distances
% A user moves 30 m. The correlation distance is 10 m, and the update
% distance is 1 m. Then, there are 3 segments (corresponding to 3 
% independent channel generations), and there are 10 snapshots in each
% segment. Thus, there are 30 snapshots in total. If the velocity is 10 m/s
% , the update time is 0.1 s, and the total moving duration is 3 s.

%% Large-scale shadow fading correlated map

% Correlation distance of SF  一个信道段的长度
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
t_dps = co_dps*t_update; % 段内t_dps总时间

%% Channel Segments
% The moving distance may not be the multiplicity of the correlation
% distance, thus the last segment may not have full length

% 一个信道段中（相关距离），大尺度衰落系数不变

% Generate the initial T-R separation distance，大小在dmin和dmax范围之间
dini = getTRSep(dmin,dmax);

% Obtain the number of channel segments 信道段的数量,
numberOfSegments = ceil(movDistance/d_co); %ceil朝正无穷向上取整

% Obtain the vector of the lengths of channel segments 
lengthOfSegments = [co_dps*ones(1,numberOfSegments-1), mod(movDistance,d_co)/d_update]; % co_dps：信道段中的采样数
if lengthOfSegments(end) == 0
    lengthOfSegments(end) = d_co/d_update;
end
% lengthOfSegments是一个向量，各元素表示各个段采样点snapshot数量
% 注意与numberOfSegments区分（--ge）

% Generate segment-wise environment parameters
% The initial T-R separation distance for each segment (the first segment 
% has been determined.)
segDist = zeros(1,numberOfSegments); 
segDist(1,1) = dini;

% The environment type and scenario type of each segment (the first segment
% has been determined.)
segEnvType = cell(1,numberOfSegments); 
segEnvType{1,1} = envType;

% The shadow fading 
segSF = zeros(1,numberOfSegments);

% The number of time clusters in each segment (independent among segments)
nTC = zeros(1,numberOfSegments); % # of time clusters in each segment

% For loop for segments

%% 关键开始
% In each loop, an anchor channel snapshot (large-scale and small-scale 
% parameters) is first generated, then, the spatial consistency procedure
% is used to update the small-scale parameters (delays, angles, powers) of
% the rest snapshots in this segment.
for segIdx = 1:numberOfSegments
    % load channel parameters
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
    SF = 1.7; mu_AOD = 1; mu_AOA = 1;X_max = 0.2; mu_tau = 123; 
    minVoidInterval = 25; sigmaCluster = 1;Gamma = 25.9; sigmaSubpath = 6; 
    gamma = 16.9; mean_ZOD = -12.6;sigma_ZOD = 5.9; std_AOD_RMSLobeAzimuthSpread = 8.5;
    std_AOD_RMSLobeElevationSpread = 2.5;distributionType_AOD = 'Gaussian'; 
    mean_ZOA = 10.8; sigma_ZOA = 5.3;std_AOA_RMSLobeAzimuthSpread = 10.5;
    std_AOA_RMSLobeElevationSpread = 11.5;distributionType_AOA = 'Laplacian';
    % RMa NLOS
    elseif strcmp(sceType,'RMa') == true && strcmp(envType,'NLOS') == true
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
    rho_mn = getIntraClusterDelays(numberOfClusterSubPaths,X_max,sceType); % rho_mn簇内子路径时延的集合（所有簇的）
    phases_mn = getSubpathPhases(rho_mn);
    tau_n = getClusterExcessTimeDelays(mu_tau,rho_mn,minVoidInterval); % 簇时延
    
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
        segDist = sqrt(sum((segInitPos-TX_Pos).^2,1)); % 某段中各路径的收发间距
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
    % Obtain shadow fading from the map. round():四舍五入
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
                                getNewPowerSpectrum(powerSpectrumOld,RFBW); % minimum time interval between adjacent subpaths in ns
    % LOS alignment （LOS修正）                       
    powerSpectrum = getLosAligned(envType,powerSpectrum);
    powerSpectrumOld = getLosAligned(envType,powerSpectrumOld);      
    
    % Generate a struct foe each independently generated channel segment
    CIR.pathDelays = powerSpectrumOld(:,1);
    pathPower = powerSpectrumOld(:,2);
    clear indNaN; indNaN = find(pathPower<=10^(Th/10));% 10^(Th/10):将dB转化为实际大小功率dBm
    pathPower(indNaN,:) = 10^(Th/10); % 将路径功率小于Th dB的置为Th
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
    CIR.environment = envType; % envType对应LOS/NLOS,会变的
    CIR.scenario = sceType; % 场景是不变的
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
    % number of multipath components in spatial and temporal domain are
    % identical. 不同
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
    %使用基于几何的方法，更新角度信息，使用多个反射面
    if (strcmp(segEnvType(segIdx),'LOS'))
        no_mpc = size(sc_powerSpectrum,1)-1;% 多径分量数量，去掉LOS
        
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
            
            % delta is the difference term between two snapshots 更新角度
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
            
            % Update power，sfMap
            sfc = sfMap(round(area/2+track(1,t+(segIdx-1)*co_dps)),round(area/2+track(2,t+(segIdx-1)*co_dps)));
            dc = sqrt(sum(track(:,t+(segIdx-1)*co_dps).^2));
            [~, Prc,~,~] = getPowerInfo(sceType,envType,f,n,sfc,TXPower,...
                                    dc,d0,p,c,u,t,RR,Pol,Fol,h_BS,folAtt,dFol);
            sc_power_cluster = getClusterPowers(sc_tau_n,Prc,Gamma,sigmaCluster,Th);
            power_temp = getSubpathPowers(sc_rho_mn,sc_power_cluster,gamma,sigmaSubpath,segEnvType(segIdx),Th);
            sc_power(:,t) = structToList(power_temp,nSP); %structToList得到各簇的子路径功率，其大小为各簇对应子路径之和。
            
            % Update phase
            sc_delay(:,t) = [sc_delay_los(1,t);sc_delay_nlos(:,t)];
            dt_delay = sc_delay(:,t)-sc_delay(:,t-1);
            sc_phase(:,t) = mod(sc_phase(:,t-1) + dt_delay*2*pi*f*1e-3, 2*pi); % 相位与时延和f相关。
            
            % Update NLOS components one by one
            for i = 1:no_mpc
                
                % Using multiplication of rotation matrices with respect to
                % x, y, z axis. 
%                 Rz1 = [cos(gcs_AOD_nlos(i,t-1)+pi), -sin(gcs_AOD_nlos(i,t-1)+pi), 0;
%                        sin(gcs_AOD_nlos(i,t-1)+pi), cos(gcs_AOD_nlos(i,t-1)+pi),0;
%                        0, 0, 1];
%                 Ry1 = [cos(pi/2-gcs_ZOD_nlos(i,t-1)), 0, sin(pi/2-gcs_ZOD_nlos(i,t-1));
%                        0, 1, 0;
%                        -sin(pi/2-gcs_ZOD_nlos(i,t-1)), 0, cos(pi/2-gcs_ZOD_nlos(i,t-1))];
%                 Ry2 = [cos(pi/2-gcs_ZOA_nlos(i,t-1)), 0, sin(pi/2-gcs_ZOA_nlos(i,t-1));
%                        0, 1, 0;
%                        -sin(pi/2-gcs_ZOA_nlos(i,t-1)), 0, cos(pi/2-gcs_ZOA_nlos(i,t-1))];
%                 Rz2 = [cos(-gcs_AOA_nlos(i,t-1)), -sin(-gcs_AOA_nlos(i,t-1)), 0;
%                        sin(-gcs_AOA_nlos(i,t-1)), cos(-gcs_AOA_nlos(i,t-1)), 0;
%                        0, 0, 1];
%                    
%                 rb = randi(2)*2-3;
%                 Rb = [1 0 0;0 rb 0;0 0 rb];
%                 R = Rz1*Ry1*Rb*Ry2*Rz2;
%                 v_RS = R*v(:,t-1);    
                
                tempBern = xBern(i);
                deltaRS = gcs_AOA_nlos(i,t-1)+(-1)^tempBern*gcs_AOD_nlos(i,t-1)+tempBern*pi;
                v_RS = mod(deltaRS+(-1)^tempBern*v(1:2,t-1),2*pi);
                v_RS = [v_RS;0];
                
                deltaAOD = v_RS'*[-sin(gcs_AOD_nlos(i,t-1));cos(gcs_AOD_nlos(i,t-1));0]*t_update;
                gcs_AOD_nlos(i,t) = gcs_AOD_nlos(i,t-1) + deltaAOD/(c*sc_delay_nlos(i,t-1)*1e-9*sin(gcs_ZOD_nlos(i,t-1)));

                deltaZOD = v_RS'*[cos(gcs_ZOD_nlos(i,t-1))*cos(gcs_AOD_nlos(i,t-1));cos(gcs_ZOD_nlos(i,t-1))*sin(gcs_AOD_nlos(i,t-1));-sin(gcs_ZOD_nlos(i,t-1))]*t_update;
                gcs_ZOD_nlos(i,t) = gcs_ZOD_nlos(i,t-1) + deltaZOD/(c*sc_delay_nlos(i,t-1)*1e-9);

                deltaAOA = v_RS'*[-sin(gcs_AOA_nlos(i,t-1));cos(gcs_AOA_nlos(i,t-1));0]*t_update;
                gcs_AOA_nlos(i,t) = gcs_AOA_nlos(i,t-1) - deltaAOA/(c*sc_delay_nlos(i,t-1)*1e-9*sin(gcs_ZOA_nlos(i,t-1)));

                deltaZOA = v_RS'*[cos(gcs_ZOA_nlos(i,t-1))*cos(gcs_AOA_nlos(i,t-1));cos(gcs_ZOA_nlos(i,t-1))*sin(gcs_AOA_nlos(i,t-1));-sin(gcs_ZOA_nlos(i,t-1))]*t_update;
                gcs_ZOA_nlos(i,t) = gcs_ZOA_nlos(i,t-1) + deltaZOA/(c*sc_delay_nlos(i,t-1)*1e-9);
                
           
            end
           
        end
        
        % Change angles back to local coordinate system
        sc_AOD_los = mod(rad2deg(pi/2 - gcs_AOD_los),360);sc_AOD_nlos = mod(rad2deg(pi/2 - gcs_AOD_nlos),360);
        sc_ZOD_los = rad2deg(pi/2 - gcs_ZOD_los);sc_ZOD_nlos = rad2deg(pi/2 - gcs_ZOD_nlos);
        sc_AOA_los = mod(rad2deg(pi/2 - gcs_AOA_los),360);sc_AOA_nlos = mod(rad2deg(pi/2 - gcs_AOA_nlos),360);
        sc_ZOA_los = rad2deg(pi/2 - gcs_ZOA_los);sc_ZOA_nlos = rad2deg(pi/2 - gcs_ZOA_nlos);
        
        % Save updated info of all snapshots
        evoCIR.pathDelays = [sc_delay_los;sc_delay_nlos];
        sc_pathPower = sc_power; % 
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
        no_mpc = size(sc_powerSpectrum,1); % NLOS场景
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
            for i = 1:no_mpc
%                 Rz1 = [cos(gcs_AOD_nlos(i,t-1)+pi), -sin(gcs_AOD_nlos(i,t-1)+pi), 0;
%                        sin(gcs_AOD_nlos(i,t-1)+pi), cos(gcs_AOD_nlos(i,t-1)+pi),0;
%                        0, 0, 1];
%                 Ry1 = [cos(pi/2-gcs_ZOD_nlos(i,t-1)), 0, sin(pi/2-gcs_ZOD_nlos(i,t-1));
%                        0, 1, 0;
%                        -sin(pi/2-gcs_ZOD_nlos(i,t-1)), 0, cos(pi/2-gcs_ZOD_nlos(i,t-1))];
%                 Ry2 = [cos(pi/2-gcs_ZOA_nlos(i,t-1)), 0, sin(pi/2-gcs_ZOA_nlos(i,t-1));
%                        0, 1, 0;
%                        -sin(pi/2-gcs_ZOA_nlos(i,t-1)), 0, cos(pi/2-gcs_ZOA_nlos(i,t-1))];
%                 Rz2 = [cos(-gcs_AOA_nlos(i,t-1)), -sin(-gcs_AOA_nlos(i,t-1)), 0;
%                        sin(-gcs_AOA_nlos(i,t-1)), cos(-gcs_AOA_nlos(i,t-1)), 0;
%                        0, 0, 1];
%                 rb = randi(2)*2-3;
%                 Rb = [1 0 0;0 rb 0;0 0 rb];
%                 R = Rz1*Ry1*Rb*Ry2*Rz2;
%                 v_RS = R*v(:,t-1);

                tempBern = xBern(i);
                deltaRS = gcs_AOA_nlos(i,t-1)+(-1)^tempBern*gcs_AOD_nlos(i,t-1)+tempBern*pi;
                v_RS = mod(deltaRS+(-1)^tempBern*v(1:2,t-1),2*pi);
                v_RS = [v_RS;0];
                
                deltaAOD = v_RS'*[-sin(gcs_AOD_nlos(i,t-1));cos(gcs_AOD_nlos(i,t-1));0]*t_update;
                gcs_AOD_nlos(i,t) = gcs_AOD_nlos(i,t-1) + deltaAOD/(c*sc_delay_nlos(i,t-1)*1e-9*sin(gcs_ZOD_nlos(i,t-1)));

                deltaZOD = v_RS'*[cos(gcs_ZOD_nlos(i,t-1))*cos(gcs_AOD_nlos(i,t-1));cos(gcs_ZOD_nlos(i,t-1))*sin(gcs_AOD_nlos(i,t-1));-sin(gcs_ZOD_nlos(i,t-1))]*t_update;
                gcs_ZOD_nlos(i,t) = gcs_ZOD_nlos(i,t-1) + deltaZOD/(c*sc_delay_nlos(i,t-1)*1e-9);

                deltaAOA = v_RS'*[-sin(gcs_AOA_nlos(i,t-1));cos(gcs_AOA_nlos(i,t-1));0]*t_update;
                gcs_AOA_nlos(i,t) = gcs_AOA_nlos(i,t-1) - deltaAOA/(c*sc_delay_nlos(i,t-1)*1e-9*sin(gcs_ZOA_nlos(i,t-1)));

                deltaZOA = v_RS'*[cos(gcs_ZOA_nlos(i,t-1))*cos(gcs_AOA_nlos(i,t-1));cos(gcs_ZOA_nlos(i,t-1))*sin(gcs_AOA_nlos(i,t-1));-sin(gcs_ZOA_nlos(i,t-1))]*t_update;
                gcs_ZOA_nlos(i,t) = gcs_ZOA_nlos(i,t-1) + deltaZOA/(c*sc_delay_nlos(i,t-1)*1e-9);

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
% Note that usually, the cluster birth and death starts from the weakest one. 
% Realize smooth transitions between two consecutive channel segments 
CIR_SISO_EVO = getTransitions(nTC,numberOfSegments,co_dps,CIR_SISO_Struct,CIR_SISO_EVO);

%% for compute
% TR3D = sqrt(sum(track.^2));
% TR2D = sqrt(sum(track(1:2,:).^2));
% omniPL = zeros(numberOfSnapshot,1);
% omniPr = zeros(numberOfSnapshot,1);
% omniDS = zeros(numberOfSnapshot,1);
% KFactor = zeros(numberOfSnapshot,1);
% 
% for i = 1:numberOfSnapshot
%     
%     % Omnidirectional channels
%     sdf = sfMap(round(area/2+track(1,i)),round(area/2+track(2,i)));
%     [PL_dB, Pr_dBm, FSPL, PLE] = getPowerInfo(sceType,envType,f,n,segSF(segIdx),TXPower,...
%                                     segDist(segIdx),d0,p,c,u,t,RR,Pol,Fol,h_BS,folAtt,dFol); 
%     omniPL(i,1) = PL_dB;
%     omniPr(i,1) = Pr_dBm;
%     CIR_tmp = CIR_SISO_EVO.(['Snapshot',num2str(i)]);
%     multipathArray = CIR_tmp.pathPowers;
%     Pr = 10*log10(multipathArray);
%     xmaxInd = find(Pr>Th); % Pr接收功率，Th
%     Pr = Pr(xmaxInd);
%     timeArray = CIR_tmp.pathDelays;
%     timeArray = timeArray(xmaxInd);
%     multipathArray = multipathArray(xmaxInd);
%     meanTau = sum(timeArray.*multipathArray)/sum(multipathArray);
%     meanTau_Sq = sum(timeArray.^2.*multipathArray)/sum(multipathArray);
%     RMSDelaySpread = sqrt(meanTau_Sq-meanTau^2);
%     omniDS(i) = RMSDelaySpread;
%     KFactor(i) = 10*log10(max(multipathArray)/(sum(multipathArray)-max(multipathArray)));
%     
% end

%% MIMO CIR for each channel snapshot
for i = 1:numberOfSnapshot
    CIR_tmp = CIR_SISO_EVO.(['Snapshot',num2str(i)]);
    [CIR_MIMO_tmp,~,~,~,~] = getLocalCIR(CIR_tmp,...
        TxArrayType,RxArrayType,Nt,Nr,Wt,Wr,dTxAnt,dRxAnt,RFBW);
    % MIMO CIR is stored
    CIR_MIMO_EVO.(['Snapshot',num2str(i)]) = CIR_MIMO_tmp;
end


toc;

