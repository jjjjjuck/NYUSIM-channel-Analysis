%% 获取信道快照

cd F:\Mat\相关矩阵距离\cmd_data
clear; 

% load all_snapshots_01
load all_snapshots_001
load all_snapshots_0001
load all_snapshots_00001
cd ..

% c1 = all_snapshots_01;
c2 = all_snapshots_001;
c3 = all_snapshots_0001;
cc = all_snapshots_00001;
data = {c2,c3,cc};
% res_cmd = 
%% 计算一段时间上的一系列CDM
for c = 1:3
    
    A = data{c}; 
    len = size(A,1)*size(A,2);
    B = reshape(A,len,size(A,3)); % 向量化,size(A,3)是时刻的样本数
    R = zeros(len,len,size(A,3));
    for i = 1:size(B,2)
        vec_H = B(:,i);
        R(:,:,i) =  vec_H*vec_H'; % 计算所有时刻的全相关矩阵
    end

    C = R(:,:,1); % 时刻1的全相关矩阵
    for i = 1:size(B,2)-1
        fenzi = trace(C*R(:,:,i+1));
        fenmu = norm(C,'fro')*norm(R(:,:,i+1),'fro'); % L2矩阵范数
        cmd(i) = 1 - fenzi./fenmu; % 时刻1 与其后所有时刻的 cmd
    end
    res_cmd(c,:)=cmd;
end


%% plot
figure1 = figure;
x = 1:size(A,3)-1;
y1 =  abs(res_cmd(1,:)); 
y2 = abs(res_cmd(2,:)); 
y3 = abs(res_cmd(3,:)); 

yy1 = smooth(y1);
yy2 = smooth(y2);
yy3 = smooth(y3);
plot(x,yy1);

hold on
plot(x,yy2);
plot(x,yy3);
xlabel('time snapshots')
ylabel('correlation matrix distance')
title('cmd between snapshots over time')
% 创建 textbox

legend1 = legend('d=0.01,movdistance=50m','d=0.001,movdistance=5m','d=0.0001,movdistance=0.5m');
ylim([0,1.1]);

set(legend1,...
    'Position',[0.547023803631465 0.226746029929509 0.338214291606631 0.0869047637212844]);

hold off

