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

function h = plotUserTrack(area,track,sfMap,sceType,envType,f,SF,dini,movDistance,velocity,FigVisibility)

%% Plot user track

figure('visible',FigVisibility);
plot(track(1,:),track(2,:),'g','LineWidth',3);
text(track(1,1)+5,track(2,1)+5,'UT','Color','m',...
    'fontsize',15,'fontweight','bold');
x1 = xlim;
y1 = ylim;
imagesc(-area/2:area/2,-area/2:area/2,sfMap);
colorbar('northoutside');
hold on;
plot(track(1,:),track(2,:),'w','LineWidth',3);
text(track(1,1)+2,track(2,1)+2,'UT','Color','w',...
    'fontsize',15,'fontweight','bold');
xlim(x1);
ylim(y1);
arrow3([track(1,1) track(2,1)],[track(1,5),track(2,5)],'LineWidth',1.5);
plot(track(1,1),track(2,1),'o','MarkerSize',10,'MarkerEdgeColor','m',...
    'MarkerFaceColor','m');
plot(0,0,'p','MarkerSize',20,'MarkerEdgeColor','y','MarkerFaceColor','y');
xlabel('Location relative to BS in X (m)','fontsize',16,'fontweight','bold');
ylabel('Location relative to BS in Y (m)','fontsize',16,'fontweight','bold');
titleName = 'User Terminal Track in Shadow Fading Map';
title(titleName,'fontsize',16,'fontweight','bold');
grid on;
grid minor;
set(gca,'YDir','normal');

yarray = linspace(0.9,0.1,8);
text_pos = 0.7;
text(text_pos,yarray(2),[num2str(f),' GHz ',char(sceType),' ',char(envType)],...
    'Units','normalized','fontsize',15,'fontweight','bold')
text(text_pos,yarray(3),[num2str(dini,'%.1f'),' m T-R Separation'],'Units','normalized',...
    'fontsize',15,'fontweight','bold')        
text(text_pos,yarray(4),['\sigma_{SF} = ', num2str(SF,'%.1f'),' dB'],...
    'Units','normalized','FontSize',15,'fontweight','bold')
text(text_pos,yarray(5),['Moving distance: ',num2str(movDistance),' m'],...
    'Units','normalized','FontSize',15,'fontweight','bold')
text(text_pos,yarray(6),['Velocity: ',num2str(velocity), ' m/s'],'Units','normalized',...
    'FontSize',15,'fontweight','bold')

h = gcf;