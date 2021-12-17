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


function h1 = plotSFMap(area,sfMap,track,sceType,envType,SF,movDistance,velocity,f,dini,FigVisibility)

%% Plot SF map
figure('visible',FigVisibility);
imagesc(-area/2:area/2,-area/2:area/2,sfMap);
colorbar('northoutside');
Xmax = max(abs(track(1,:)))+10;
Ymax = max(abs(track(2,:)))+10;
hold on;
plot(track(1,:),track(2,:),'w','LineWidth',3);
titleName = 'Intensity Map of Spatially Correlated Shadow Fading (dB)';

title(titleName,'fontsize',16,'fontweight','bold');
% hold on;
plot(track(1,1),track(2,1),'o','MarkerSize',5,'MarkerEdgeColor','m',...
    'MarkerFaceColor','m');
plot(0,0,'p','MarkerSize',20,'MarkerEdgeColor','y','MarkerFaceColor','y');
% plot([0 track(1,1)],[0 track(2,1)],':','LineWidth',2);
grid on;
grid minor;

xlabel('X (m)','fontsize',16,'fontweight','bold');
ylabel('Y (m)','fontsize',16,'fontweight','bold');
text(8,8,'BS','Color','y',...
    'fontsize',15,'fontweight','bold');

text(track(1,1)+5,track(2,1)+5,'UT','Color','m',...
    'fontsize',15,'fontweight','bold');
xlim([-Xmax,Xmax]);
ylim([-Ymax,Ymax]);
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
% text(text_pos,yarray(7),['PLE = ', num2str(PLE,'%.1f')],'Units','normalized',...
%     'FontSize',15,'fontweight','bold')

h1 = gcf;

