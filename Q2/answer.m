%% fit验证
clear;clc;
load('附件2-风电机组采集数据.mat','data_TS_WF');
WT1 = data_TS_WF.WF_1.WT; 
datat = 1901;
WT=WT1;
%
Tshaft_data = zeros(100,100);
Tshaft_m = zeros(2,100);
Ft_data = zeros(100,100);
Ft_m = zeros(2,100);
for wind_num = 1:100
    tempwt = WT{1,wind_num};
    Pref = tempwt.inputs(:,1);
    Wind = tempwt.inputs(:,2);
    Tshaft = tempwt.outputs(:,1);
    Ft = tempwt.outputs(:,2);
    Pitch = tempwt.states(:,1);
    Wr = tempwt.states(:,2);
    %
    X1 = Pref(datat-1:end-1);
    Y1 = Tshaft(datat:end);
    % 假设 Y 是列向量, X1-X6 是特征矩阵中的列
    Y2 = [Ft(datat:end)];  % 你的目标变量Y
    X2 = [Wind(datat:end) Pitch(datat:end) Wr(datat-1:end-1)];  % 特征变量矩阵
    
    load fitTshaft
    Y1pre = fitTshaft(2)+fitTshaft(3).*X1;
    Tshaft_data(:,wind_num) = Y1pre;
    
    real = Y1;
    pre = Y1pre;
    mse = sum((real - pre).^2) / length(real);
    mape = sum(abs((real - pre) ./ real)) / length(real) * 100;
    Tshaft_m(1,wind_num) = mse;
    Tshaft_m(2,wind_num) = mape;
    % figure()
    % plot(1901:2000,Y1)
    % hold on
    % plot(1901:2000,Y1pre)
    % legend('real','pre')
    % title('主轴转矩估计')
    load fitFt
    Y2pre = fitFt(2)+fitFt(3).*X2(:,1)+fitFt(4).*X2(:,2)+fitFt(5).*X2(:,3);
    Ft_data(:,wind_num) = Y2pre;
    real = Y2;
    pre = Y2pre;
    mse = sum((real - pre).^2) / length(real);
    mape = sum(abs((real - pre) ./ real)) / length(real) * 100;
    Ft_m(1,wind_num) = mse;
    Ft_m(2,wind_num) = mape;
% figure()
% plot(1901:2000,Y2)
% hold on
% plot(1901:2000,Y2pre)
% legend('real','pre')
% title('推力估计')
end
%%
figure()
plot(1:100,Tshaft_m(1,:))
hold on 
title('Tshaft mse')
figure()
plot(1:100,Tshaft_m(2,:))
hold on 
title('Tshaft mape')
%%
figure()
plot(1:100,Ft_m(1,:))
hold on 
title('Ft mse')
figure()
plot(1:100,Ft_m(2,:))
hold on 
title('Ft mape')
%%
figure()
boxplot(Tshaft_m(2,:))
title('Tshaft mape')

figure()
boxplot(Ft_m(2,:))
title('Ft mape')

figure()
histogram(Tshaft_m(2,:))
title('Tshaft mape')

figure()
histogram(Ft_m(2,:))
title('Ft mape')