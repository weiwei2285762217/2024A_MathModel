%% fit
datat = 101;
X1 = [Wind(1:datat-2) Pref(1:datat-2)];  % 特征变量矩阵
Y1 = Pitch(3:datat);

X2 = [Wind(1:datat-1) Pref(1:datat-1)];
Y2 = Wr(2:datat);

X=X1;
Y=Y1;
% 在特征矩阵前添加常数项 (用于拟合常数项 β0)
X_with_intercept = [ones(size(X,1),1), X];  % 添加一列全1
% 执行线性回归拟合
mdl = fitlm(X_with_intercept, Y);
fitPitch = mdl.Coefficients.Estimate;
save('fitPitch.mat','fitPitch')
disp(mdl);
%
X=X2;
Y=Y2;
% 在特征矩阵前添加常数项 (用于拟合常数项 β0)
X_with_intercept = [ones(size(X,1),1), X];  % 添加一列全1

% 执行线性回归拟合
mdl = fitlm(X_with_intercept, Y);
fitWr = mdl.Coefficients.Estimate;
save('fitWr.mat','fitWr')
disp(mdl);
% 输出拟合结果

%% fit验证
datat = 1901;

X1 = [Wind(datat-2:end-2) Pref(datat-2:end-2)];  % 特征变量矩阵
Y1 = Pitch(datat:end);

X2 = [Wind(datat-1:end-1) Pref(datat-1:end-1)];
Y2 = Wr(datat:end);

load fitPitch
Y1pre = fitPitch(2)+fitPitch(3).*X1(:,1)+fitPitch(4).*X1(:,2);
real = Y1;
pre = Y1pre;
figure()
plot(1901:2000,Y1)
hold on
plot(1901:2000,Y1pre)
legend('real','pre')
title('Pitch估计')
mse = sum((real - pre).^2) / length(real)
mape = sum(abs((real - pre) ./ real)) / length(real) * 100

load fitWr
Y2pre = fitWr(2)+fitWr(3).*X2(:,1)+fitWr(4).*X2(:,2);
real = Y2;
pre = Y2pre;
mse = sum((real - pre).^2) / length(real)
mape = sum(abs((real - pre) ./ real)) / length(real) * 100
figure()
plot(1901:2000,Y2)
hold on
plot(1901:2000,Y2pre)
legend('real','pre')
title('Wr估计')
%% ft_test
datat = 1901;

X1 = [Wind(datat-2:end-2) Pref(datat-2:end-2)];  % 特征变量矩阵
Y1 = Pitch(datat:end);

X2 = [Wind(datat-1:end-1) Pref(datat-1:end-1)];
Y2 = Wr(datat:end);

Y3 = [Ft(datat:end)];  % 你的目标变量Y
X3 = [Wind(datat:end) Pitch(datat:end) Wr(datat-1:end-1)];  % 特征变量矩阵

load fitFt
Y3pre = fitFt(2)+fitFt(3).*X3(:,1)+fitFt(4).*X3(:,2)+fitFt(5).*X3(:,3);

load fitPitch
Y1pre = fitPitch(2)+fitPitch(3).*X1(:,1)+fitPitch(4).*X1(:,2);
load fitWr
Y2pre = fitWr(2)+fitWr(3).*X1(:,1)+fitWr(4).*X1(:,2);
Y3pre = fitFt(2)+fitFt(3).*X3(:,1)+fitFt(4).*Y1pre+fitFt(5).*Y2pre;
real = Y3;
pre = Y3pre;
mse = sum((real - pre).^2) / length(real)
mape = sum(abs((real - pre) ./ real)) / length(real) * 100
figure()
plot(1901:2000,Y3)
hold on
plot(1901:2000,Y3pre)
legend('real','pre')
title('推力估计')