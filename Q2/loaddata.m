%% loaddata
clear;clc;
load('附件2-风电机组采集数据.mat','data_TS_WF');
WT1 = data_TS_WF.WF_1.WT; 
WT2 = data_TS_WF.WF_2.WT;
%%
WT=WT1;
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
%% test
figure()
plot(Time,Ft)
hold
plot(Time,2.*Pref./Wind)
legend('real','pre')
xlim([0 100])

figure()
plot(Time,Tshaft)
hold
plot(Time,30*Pref./Wr/pi)
legend('real','pre')
xlim([0 100])