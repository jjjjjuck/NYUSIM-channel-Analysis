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

function h1 = plotConDirPDP(CIR,theta_3dB_TX,phi_3dB_TX,theta_3dB_RX,phi_3dB_RX,Th,f,sceType,envType,dini,d_update,d_co,movDistance,velocity,FigVisibility)

h1 = figure('visible',FigVisibility);
FS = 16;
no_snap = numel(fieldnames(CIR));
maxDelay = 0;
minDelay = 1000;
maxPower = Th;

for i = 1:no_snap
    CIR_tmp = CIR.(['Snapshot',num2str(i)]);
    ps = [CIR_tmp.pathDelays, CIR_tmp.pathPowers, CIR_tmp.pathPhases,...
    CIR_tmp.AODs, CIR_tmp.ZODs, CIR_tmp.AOAs, CIR_tmp.ZOAs];

    theta_TX_d = ps(:,4);
    phi_TX_d = ps(:,5);
    theta_RX_d = ps(:,6);
    phi_RX_d = ps(:,7);
    nPath = size(ps,1);

    [~, maxIndex] = max(ps(:,2));

    [TX_Dir_Gain_Mat, RX_Dir_Gain_Mat, G_TX, G_RX] = getDirectiveGains(theta_3dB_TX,...
        phi_3dB_TX,theta_3dB_RX,phi_3dB_RX,theta_TX_d(maxIndex),phi_TX_d(maxIndex),...
        theta_RX_d(maxIndex),phi_RX_d(maxIndex),ps);

    [timeArray_Dir, multipathArray_Dir] = getDirPDP(ps,...
        TX_Dir_Gain_Mat,RX_Dir_Gain_Mat);
    
    Pr = 10*log10(multipathArray_Dir);
    xmaxInd = find(Pr>Th);
    Pr = Pr(xmaxInd);
    timeArray_Dir = timeArray_Dir(xmaxInd);
    maxDelayTemp = max(timeArray_Dir);
    minDelayTemp = min(timeArray_Dir);
    maxPowerTemp = max(Pr);
        
    if maxDelayTemp > maxDelay
        maxDelay = maxDelayTemp;
    end
    if minDelayTemp < minDelay
        minDelay = minDelayTemp;
    end
    if maxPowerTemp > maxPower
        maxPower = maxPowerTemp;
    end
    
    if numel(xmaxInd) == 1
        RMSDelaySpread=0;
    else
    end
    
    if ~isempty(xmaxInd)
        stem3(ones(length(timeArray_Dir),1)*i,timeArray_Dir,Pr,'BaseValue',Th,'LineStyle','-','Marker','none','linewidth',1);
        hold on;
        set(gca,'Ydir','reverse'); 
    else
        text(0.18,0.6,['No Detectable Multipath Components',char(10),'above the Threshold of ',num2str(Th),' dBm'],...
                'Units','normalized','fontsize',15,'fontweight','bold');
    end
    
    
end
titleName = 'Directional Power Delay Profile (PDP) Evolution with Strongest Power';
title(titleName,'fontsize',20,'fontweight','bold')

xlabel('Track Distance (m)','fontsize',FS); 
ylabel('Absolute Propagation Time (ns)','fontsize',FS);
zlabel('Received Power (dBm)','fontsize',FS);
xh = get(gca,'XLabel'); set(xh, 'Units', 'Normalized'); 
pos = get(xh, 'Position'); set(xh, 'Position',pos.*[1.1,-0.5,0],'Rotation',15);
yh = get(gca,'YLabel'); set(yh, 'Units', 'Normalized'); 
pos = get(yh, 'Position'); set(yh, 'Position',pos.*[0.7,-0.6,0],'Rotation',-25);
set(gca,'fontsize',FS); 
xlim([1 no_snap]);
ylim([.985*minDelay 1.2*maxDelay]);
zlim([Th maxPower+5]);
xt = xticks;
xticks([1 xt]);
yarray = linspace(0.9,0.1,8);
text_pos = 0.8;
text(text_pos,yarray(2),[num2str(f),' GHz ',char(sceType),' ',char(envType)],...
    'Units','normalized','fontsize',15,'fontweight','bold')
text(text_pos,yarray(3),[num2str(dini,'%.1f'),' m T-R Separation'],'Units','normalized',...
    'fontsize',15,'fontweight','bold')        
text(text_pos,yarray(4),['Update distance = ', num2str(d_update,'%.1f'),' m'],...
    'Units','normalized','FontSize',15,'fontweight','bold')
text(text_pos,yarray(5),['Track distance: ',num2str(movDistance,'%.1f'),' m'],...
    'Units','normalized','FontSize',15,'fontweight','bold')
text(text_pos,yarray(6),['Velocity: ',num2str(velocity), ' m/s'],'Units','normalized',...
    'FontSize',15,'fontweight','bold')
text(text_pos,yarray(7),['Segment Length = ', num2str(d_co,'%.1f'),' m'],...
    'Units','normalized','FontSize',15,'fontweight','bold')
set(gca,'fontsize',FS);

end
