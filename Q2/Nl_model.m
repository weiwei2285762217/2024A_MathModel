%% loaddata
clear;clc;
load('附件2-风电机组采集数据.mat','data_TS_WF');
WT1 = data_TS_WF.WF_1.WT; 
WT2 = data_TS_WF.WF_2.WT;
%%
WT=WT1;
WTnum = 2;
Time = WT{1, WTnum}.time;
Dt = Time(2)-Time(1);
Pref = WT{1,WTnum}.inputs(:,1);
Wind = WT{1, WTnum}.inputs(:,2);
Tshaft = WT{1, WTnum}.outputs(:,1);
Ft = WT{1, WTnum}.outputs(:,2);
Pout = WT{1, WTnum}.outputs(:,3);
Pitch = WT{1, WTnum}.states(:,1);
Wr = WT{1, WTnum}.states(:,2)*pi/30;%rpm-rad/s
Wg= WT{1, WTnum}.states(:,3)*pi/30;


%% Tshaft = Ft - Tr*Wr_a
Wr_p = zeros(length(Time),1);
Wr_p(1) = Wr(1);
% Jr = 38759236e6;
Jr = 1e10;
for i=2:Time(10)
    Ta = Pout(i)./Wr(i);
    Wr_p(i) = (Dt./Jr)*(Ta-Tshaft(i))+Wr_p(i-1);
    
end

figure()
plot(Time,Wr_p)
hold on
plot(Time,Wr)
legend('pre','real')
title('低速轴转速预测')
xlim([0 10])
mse = sum((Wr - Wr_p).^2) / length(Wr)
mape = sum(abs((Wr - Wr_p) ./ Wr)) / length(Wr) * 100
%%
Jr = zeros(2000,1);
for i=2:Time(end)
    Jr(i) = Dt*(Pref(i-1)./Wr(i-1)-Tshaft(i-1))./(Wr(i)-Wr(i-1));

end
figure()
plot(Time,Jr)