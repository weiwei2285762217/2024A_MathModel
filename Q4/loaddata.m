%% loaddata
clear;clc;
load('附件3-噪声和延迟作用下的采集数据.mat')
WT2 = data_TS_WF.WF_1.WT; 
load('附件2-风电机组采集数据.mat')
WT1 = data_TS_WF.WF_1.WT;
%%
WT=WT2;
tempwt = WT{1,1};
Time = tempwt.time;
Dt = Time(2)-Time(1);
Pref = tempwt.inputs(:,1);
Wind = tempwt.inputs(:,2);
Tshaft = tempwt.outputs(:,1);
Ft = tempwt.outputs(:,2);
Pout = tempwt.outputs(:,3);
Pitch = tempwt.states(:,1);
Wr = tempwt.states(:,2);
Wg= tempwt.states(:,3);
WTdata = [Pref Wind Tshaft Ft Pout Pitch Wr Wg];
%% view
figure()
plot(Time,Pref)
hold on
plot(Time,Pout)
title('功率')
legend('Pref','Pout')
xlim([0 100])
%% 功率
figure()
yyaxis left
plot(Time,Pref)
hold on
title('功率-主轴')

yyaxis right
plot(Time(1:end-1),Tshaft(2:end))

legend('Pref','Tshaft')
xlim([0 100])
%% 风速-桨距角
figure()
yyaxis left
plot(Time,Wind)
hold on
title('风速-桨距角')


yyaxis right
plot(Time(1:end-1),Pitch(2:end))

legend('wind','pitch')
xlim([0 100])

%% 风速-推力
figure()
yyaxis left
plot(Time,Wind)
hold on
plot(Time,Pitch)
title('风速-推力-桨距角')


yyaxis right
plot(Time,Ft)

legend('wind','pitch','Ft')
xlim([0 100])
%% 推力-桨距角
figure()
yyaxis left
plot(Time,Ft)
hold on
title('推力-桨距角')


yyaxis right
plot(Time,Pitch)

legend('Ft','pitch')
xlim([0 100])

%% 转速
figure()
N=97;
plot(Time,Wr*N)
hold on
plot(Time,Wg)
legend('Wr','Wg')
title('转速')
xlim([0 100])
%% corr1
% 示例数据
data = WTdata; % 10行5列的随机数据，替换为你的实际数据

% 计算相关性矩阵
correlationMatrix = corr(data);
%Pref Wind Tshaft Ft Pout Pitch Wr Wg
xLabels = {'给定功率', '风速', '主轴转矩', '推力', '输出功率', '桨距角','低速转速','高速转速'};
yLabels = {'给定功率', '风速', '主轴转矩', '推力', '输出功率', '桨距角','低速转速','高速转速'};
% 创建热力图
figure;
heatmap(xLabels, yLabels, correlationMatrix, 'ColorLimits', [-1 1], 'Colormap', jet);

% 设置图形的标签
xlabel('变量');
ylabel('变量');
title('相关性热力图');
%% corr-delay
% 示例数据
% WTdata(:,2) = WTdata(:,2).^2;
data = [WTdata(1:end-1,:) WTdata(2:end,:)]; % 10行5列的随机数据，替换为你的实际数据

% 计算相关性矩阵
correlationMatrix = corr(data);
%Pref Wind Tshaft Ft Pout Pitch Wr Wg
xLabels = {'给定功率', '风速', '主轴转矩', '推力', '输出功率', '桨距角','低速转速','高速转速',...
    '给定功率2', '风速2', '主轴转矩2', '推力2', '输出功率2', '桨距角2','低速转速2','高速转速2'};
yLabels = {'给定功率', '风速', '主轴转矩', '推力', '输出功率', '桨距角','低速转速','高速转速',...
    '给定功率2', '风速2', '主轴转矩2', '推力2', '输出功率2', '桨距角2','低速转速2','高速转速2'};
% 创建热力图
figure;
heatmap(xLabels, yLabels, correlationMatrix, 'ColorLimits', [-1 1], 'Colormap', jet);

% 设置图形的标签
xlabel('变量');
ylabel('变量');
title('相关性热力图');
%% PCA
datat = 101;
% 假设你的数据矩阵为 data，行是样本，列是变量
Y2 = [Ft(2:datat)];  % 你的目标变量Y
X2 = [Wind(2:datat) Pitch(2:datat) Wr(1:datat-1)];  % 特征变量矩阵
X2_1 = [Wind(2:datat) Pitch(2:datat) Wr(1:datat-1) Tshaft(1:datat-1) Ft(1:datat-1)];  
X2_2 = [Wind(2:datat) Pitch(2:datat) Wr(1:datat-1) Tshaft(1:datat-1)]; 
X2_3 = [Wind(2:datat) Pitch(2:datat) Wr(1:datat-1) Ft(1:datat-1)]; 
data = [X2];
% 步骤1：标准化数据
data_standardized = zscore(data);

% 步骤2：进行 PCA
[coeff, score, latent, tsquared, explained, mu] = pca(data_standardized);

% coeff：主成分系数
% score：主成分得分
% latent：每个主成分的特征值
% explained：每个主成分解释的方差百分比

% 步骤3：可视化前两个主成分
figure;
scatter(score(:, 1), score(:, 2));
xlabel('主成分 1');
ylabel('主成分 2');
title('PCA 主成分得分图');

% 可选：显示每个主成分解释的方差百分比
figure;
pareto(explained);
xlabel('主成分');
ylabel('解释的方差百分比');
title('主成分方差贡献');
%% fit
datat = 101;
X1 = Pref(1:datat-1);
Y1 = Tshaft(2:datat);
% 假设 Y 是列向量, X1-X6 是特征矩阵中的列
Y2 = [Ft(2:datat)];  % 你的目标变量Y
X2 = [Wind(2:datat) Pitch(2:datat) Wr(1:datat-1)];  % 特征变量矩阵
X2_1 = [Wind(2:datat) Pitch(2:datat) Wr(1:datat-1) Tshaft(1:datat-1) Ft(1:datat-1)];  % 特征变量矩阵

X=X1;
Y=Y1;
% 在特征矩阵前添加常数项 (用于拟合常数项 β0)
X_with_intercept = [ones(size(X,1),1), X];  % 添加一列全1

% 执行线性回归拟合
mdl = fitlm(X_with_intercept, Y);
fitTshaft = mdl.Coefficients.Estimate;
save('fitTshaft.mat','fitTshaft')
disp(mdl);
X=X2;
Y=Y2;
% 在特征矩阵前添加常数项 (用于拟合常数项 β0)
X_with_intercept = [ones(size(X,1),1), X];  % 添加一列全1

% 执行线性回归拟合
mdl = fitlm(X_with_intercept, Y);
fitFt = mdl.Coefficients.Estimate;
save('fitFt.mat','fitFt')
disp(mdl);
% 输出拟合结果
% disp(mdl);

% % 显示回归系数
% coefficients = mdl.Coefficients.Estimate; % 获取回归系数
% beta_0 = coefficients(1);  % 常数项
% beta_1 = coefficients(2);  % X1 的系数
% beta_2 = coefficients(3);  % X2 的系数
% beta_3 = coefficients(4);  % X3 的系数
% beta_4 = coefficients(5);  % X4 的系数
% beta_5 = coefficients(6);  % X5 的系数
% beta_6 = coefficients(7);  % X6 的系数
% 
% % 显示回归方程
% fprintf('Y = %.4f + %.4f * X1 + %.4f * X2 + %.4f * X3 + %.4f * X4 + %.4f * X5 + %.4f * X6\n', ...
%     beta_0, beta_1, beta_2, beta_3, beta_4, beta_5, beta_6);
%% fit验证
datat = 1901;
X1 = Pref(datat-1:end-1);
Y1 = Tshaft(datat:end);
% 假设 Y 是列向量, X1-X6 是特征矩阵中的列
Y2 = [Ft(datat:end)];  % 你的目标变量Y
X2 = [Wind(datat:end) Pitch(datat:end) Wr(datat-1:end-1)];  % 特征变量矩阵
X2_1 = [Wind(datat:end) Pitch(datat:end) Wr(datat-1:end-1) Tshaft(datat-1:end-1) Ft(datat-1:end-1)];  % 特征变量矩阵

load fitTshaft
Y1pre = fitTshaft(2)+fitTshaft(3).*X1;
real = Y1;
pre = Y1pre;
figure()
plot(1901:2000,Y1)
hold on
plot(1901:2000,Y1pre)
legend('real','pre')
title('主轴转矩估计')
mse = sum((real - pre).^2) / length(real)
mape = sum(abs((real - pre) ./ real)) / length(real) * 100
load fitFt
Y2pre = fitFt(2)+fitFt(3).*X2(:,1)+fitFt(4).*X2(:,2)+fitFt(5).*X2(:,3);
real = Y2;
pre = Y2pre;
mse = sum((real - pre).^2) / length(real)
mape = sum(abs((real - pre) ./ real)) / length(real) * 100
figure()
plot(1901:2000,Y2)
hold on
plot(1901:2000,Y2pre)
legend('real','pre')
title('推力估计')
%% corr-^2
% 示例数据
% % WTdata(:,2) = WTdata(:,2).^2;
% data = [WTdata(1:end-1,:) log(WTdata(1:end-1,:))]; % 10行5列的随机数据，替换为你的实际数据
% 
% % 计算相关性矩阵
% correlationMatrix = corr(data);
% %Pref Wind Tshaft Ft Pout Pitch Wr Wg
% xLabels = {'给定功率', '风速', '主轴转矩', '推力', '输出功率', '桨距角','低速转速','高速转速',...
%     '给定功率2', '风速2', '主轴转矩2', '推力2', '输出功率2', '桨距角2','低速转速2','高速转速2'};
% yLabels = {'给定功率', '风速', '主轴转矩', '推力', '输出功率', '桨距角','低速转速','高速转速',...
%     '给定功率2', '风速2', '主轴转矩2', '推力2', '输出功率2', '桨距角2','低速转速2','高速转速2'};
% % 创建热力图
% figure;
% heatmap(xLabels, yLabels, correlationMatrix, 'ColorLimits', [-1 1], 'Colormap', jet);
% 
% % 设置图形的标签
% xlabel('变量');
% ylabel('变量');
% title('相关性热力图');

