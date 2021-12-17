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

function h1 = plotConPDP(CIR,Th,f,sceType,envType,dini,d_update,d_co,movDistance,velocity,FigVisibility)

h1 = figure('visible',FigVisibility);
FS = 16;
no_snap = numel(fieldnames(CIR));
maxDelay = 0;
minDelay = 1000;
maxPower = Th;

for i = 1:no_snap
    CIR_tmp = CIR.(['Snapshot',num2str(i)]);
    multipathArray = CIR_tmp.pathPowers;
    Pr = 10*log10(multipathArray);
    xmaxInd = find(Pr>Th);
    Pr = Pr(xmaxInd);
    timeArray = CIR_tmp.pathDelays;
    timeArray = timeArray(xmaxInd);
    maxDelayTemp = max(timeArray);
    minDelayTemp = min(timeArray);
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
%     Pr_Lin = sum(multipathArray);
%     meanTau = sum(timeArray.*multipathArray)/sum(multipathArray);
%     meanTau_Sq=sum(timeArray.^2.*multipathArray)/sum(multipathArray);
%     RMSDelaySpread = sqrt(meanTau_Sq-meanTau^2);
%     
    if numel(xmaxInd) == 1
        RMSDelaySpread=0;
    else
    end

    if ~isempty(xmaxInd)
        stem3(ones(length(timeArray),1)*i,timeArray,Pr,'BaseValue',Th,'LineStyle','-','Marker','none','linewidth',1);
        hold on;
        set(gca,'Ydir','reverse'); 
        
    else
            text(0.18,0.6,['No Detectable Multipath Components',char(10),'above the Threshold of ',num2str(Th),' dBm'],...
                'Units','normalized','fontsize',15,'fontweight','bold');
    end
end
titleName = 'Omnidirectional Power Delay Profile (PDP) Evolution';

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
ylim([.985*minDelay 1.4*maxDelay]);
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