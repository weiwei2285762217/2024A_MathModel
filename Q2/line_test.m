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
data = [WTdata(1:end-3,:) WTdata(4:end,:)]; % 10行5列的随机数据，替换为你的实际数据

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
X1_1 = [Pref(1:datat-1) ];
Y1 = Tshaft(2:datat);
% 假设 Y 是列向量, X1-X6 是特征矩阵中的列
Y2 = [Ft(2:datat)];  % 你的目标变量Y
X2 = [Wind(2:datat) Pitch(2:datat) Wr(1:datat-1)];  % 特征变量矩阵
X2_1 = [Wind(2:datat) Pitch(2:datat) Wr(1:datat-1) Tshaft(1:datat-1) Ft(1:datat-1)];  % 特征变量矩阵
X3 = [Ft(1:datat-1)];
X=X1;
Y=Y1;
% 在特征矩阵前添加常数项 (用于拟合常数项 β0)
X_with_intercept = [ones(size(X,1),1), X];  % 添加一列全1
% 执行线性回归拟合
mdl = fitlm(X_with_intercept, Y);
fitTshaft = mdl.Coefficients.Estimate;
save('fitTshaft.mat','fitTshaft')
disp(mdl);
%
X=X2;
Y=Y2;
% 在特征矩阵前添加常数项 (用于拟合常数项 β0)
X_with_intercept = [ones(size(X,1),1), X];  % 添加一列全1

% 执行线性回归拟合
mdl = fitlm(X_with_intercept, Y);
fitFt = mdl.Coefficients.Estimate;
save('fitFt.mat','fitFt')
disp(mdl);
% % 输出拟合结果
% X=X2_1;
% Y=Y2;
% % 在特征矩阵前添加常数项 (用于拟合常数项 β0)
% X_with_intercept = [ones(size(X,1),1), X];  % 添加一列全1
% 
% % 执行线性回归拟合
% mdl = fitlm(X_with_intercept, Y);
% fitFt2 = mdl.Coefficients.Estimate;
% save('fitFt2.mat','fitFt2')
% disp(mdl);
% 
% X=X3;
% Y=Y2;
% % 在特征矩阵前添加常数项 (用于拟合常数项 β0)
% X_with_intercept = [ones(size(X,1),1), X];  % 添加一列全1
% 
% % 执行线性回归拟合
% mdl = fitlm(X_with_intercept, Y);
% fitFt3 = mdl.Coefficients.Estimate;
% save('fitFt3.mat','fitFt3')
% disp(mdl);
%% fit验证
datat = 1901;
X1 = Pref(datat-1:end-1);
Y1 = Tshaft(datat:end);
% 假设 Y 是列向量, X1-X6 是特征矩阵中的列
Y2 = [Ft(datat:end)];  % 你的目标变量Y
X2 = [Wind(datat:end) Pitch(datat:end) Wr(datat-1:end-1)];  % 特征变量矩阵
X2_1 = [Wind(datat:end) Pitch(datat:end) Wr(datat-1:end-1) Tshaft(datat-1:end-1) Ft(datat-1:end-1)];  % 特征变量矩阵
X3 = [Wind(datat:end) Ft(datat-1:end-1)];
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
% 
% load fitFt2
% Y2_1pre = fitFt2(1)+fitFt2(3).*X2_1(:,1)+fitFt2(4).*X2_1(:,2)+fitFt2(5).*X2_1(:,3)+fitFt2(6).*X2_1(:,4)+fitFt2(7).*X2_1(:,5);
% real = Y2;
% pre = Y2_1pre;
% mse = sum((real - pre).^2) / length(real)
% mape = sum(abs((real - pre) ./ real)) / length(real) * 100
% figure()
% plot(1901:2000,Y2)
% hold on
% plot(1901:2000,Y2_1pre)
% legend('real','pre')
% title('推力估计')
% 
% load fitFt3
% Y3pre = 1;
% real = Y2;
% pre = Y3pre;
% mse = sum((real - pre).^2) / length(real)
% mape = sum(abs((real - pre) ./ real)) / length(real) * 100
% figure()
% plot(1901:2000,Y2)
% hold on
% plot(1901:2000,Y3pre)
% legend('real','pre')
% title('推力估计')
%% corr-^2
% 示例数据
% WTdata(:,2) = WTdata(:,2).^2;
data = [WTdata(1:end-1,:) log(WTdata(1:end-1,:))]; % 10行5列的随机数据，替换为你的实际数据
% data = [WTdata(1:end-1,:) exp(WTdata(1:end-1,:))]; 
data = [WTdata(1:end-1,:) WTdata(2:end,:)]; 
data = [WTdata(1:end-2,:) WTdata(3:end,:)]; 
data = [WTdata(1:end-3,:) WTdata(4:end,:)]; 
% data = [WTdata(1:end-1,:) (WTdata(1:end-1,:)).^2];
% data = [WTdata(1:end-1,:) 1./(WTdata(1:end-1,:))];
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
