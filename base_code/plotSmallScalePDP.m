%%% NYUSIM - User License %%%

% Copyright (c) 2016-2021 New York University and NYU WIRELESS

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

function [h,X,Y,Pr_H] = plotSmallScalePDP(FigVisibility,CIR_MIMO,Nr,Th,TXPower,f,RFBW,...
    sceType,envType,TRDistance)

figure('Visible',FigVisibility);
% Fontsize
FS = 12; 
% Time delay of multipath components
timeDelay = CIR_MIMO.pathDelays; 
% PDPs along RX antenna elements
Pr_H = 10.*log10((abs(CIR_MIMO.HSmallScale)).^2) + TXPower; 
% If the received power is smaller than the threshold, set it to the threshold
Pr_H(Pr_H<Th) = Th;

X = (repmat([1:Nr],length(timeDelay),1)-1)./2; Y = repmat(timeDelay,1,Nr);
if max(max(Pr_H)) > Th
stem3(X,Y,Pr_H,'BaseValue',Th,'LineStyle','-','Marker','none',...
    'linewidth',1.5,'color','b');
else
    text(0.18,0.6,['No Detectable Multipath Components',char(10),'above the Threshold of ',num2str(Th),' dBm'],...
        'Units','normalized','fontsize',15,'fontweight','bold');
end
xlabel('Increment (# of Wavelengths)','fontsize',FS); 
ylabel('Propagation Time Delay (ns)','fontsize',FS);
zlabel('Received Power (dBm)','fontsize',FS);
xh = get(gca,'XLabel'); set(xh, 'Units', 'Normalized'); 
pos = get(xh, 'Position'); set(xh, 'Position',pos.*[1.1,-0.5,0],'Rotation',15);
yh = get(gca,'YLabel'); set(yh, 'Units', 'Normalized'); 
pos = get(yh, 'Position'); set(yh, 'Position',pos.*[0.7,-0.6,0],'Rotation',-25);
set(gca,'Ydir','reverse'); 
set(gca,'fontsize',FS); 
zlim([Th max(max(max(Pr_H)),Th)+1]);
title(['Small Scale PDPs - ',num2str(f),' GHz, ' num2str(RFBW),' MHz, ', sceType,' ',envType,' ',...
    num2str(TRDistance,'%.1f'),' m T-R Separation']); 
h = gcf;  