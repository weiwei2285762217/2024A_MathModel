%% loaddata
clear;clc;
load('附件3-噪声和延迟作用下的采集数据.mat')
WT2 = data_TS_WF.WF_1.WT; 
load('附件2-风电机组采集数据.mat')
WT1 = data_TS_WF.WF_1.WT;
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
%%
WT=WT2;
tempwt2 = WT{1,1};
Time2 = tempwt2.time;
Dt2 = Time(2)-Time(1);
Pref2 = tempwt2.inputs(:,1);
Wind2 = tempwt2.inputs(:,2);
Tshaft2 = tempwt2.outputs(:,1);
Ft2 = tempwt2.outputs(:,2);
Pout2 = tempwt2.outputs(:,3);
Pitch2 = tempwt2.states(:,1);
Wr2 = tempwt2.states(:,2);
Wg2 = tempwt2.states(:,3);
WTdata2 = [Pref2 Wind2 Tshaft2 Ft2 Pout2 Pitch2 Wr2 Wg2];
%% output
figure()
plot(Time,Ft)
hold on
plot(Time,Ft2)
title('Ft')
legend('Ft','zaosheng')
xlim([0 100])


figure()
plot(Time,Tshaft)
hold on
plot(Time,Tshaft2)
title('Tshaft')
legend('Tshaft','zaosheng')
xlim([0 100])

figure()
plot(Time,Pout)
hold on
plot(Time,Pout2)
title('功率')
legend('Pref','zaosheng')
xlim([0 100])
%% input
figure()
plot(Time,Pref)
hold on
plot(Time,Pref2)
title('功率')
legend('Pref','zaosheng')
xlim([0 100])

figure()
plot(Time,Wind)
hold on
plot(Time,Wind2)
title('风速')
legend('Wind','zaosheng')
xlim([200 300])
%找延迟
index = find(abs(Wind2(1:end-1)-Wind2(2:end))==0)

% figure()
% plot(Time(2:end),abs(Wind2(1:end-1)-Wind2(2:end)))
% title('延迟')
%% state
figure()
plot(Time,Pitch)
hold on
plot(Time,Pitch2)
title('桨距角')
legend('Pitch','zaosheng')
xlim([0 100])
figure()
plot(Time,Wg)
hold on
plot(Time,Wg2)
title('Wr')
legend('Wr','zaosheng')
xlim([0 100])
figure()
plot(Time,Wr)
hold on
plot(Time,Wr2)
title('Wr')
legend('Wr','zaosheng')
xlim([0 100])
%% 延迟分析
wind_num = 1;
for wind_num=1:100
    WT=WT2;
    tempwt2 = WT{1,wind_num};
    Pref2_all(:,wind_num) = tempwt2.inputs(:,1);
    Wind2_all(:,wind_num) = tempwt2.inputs(:,2);
    Tshaft2_all(:,wind_num) = tempwt2.outputs(:,1);
    Ft2_all(:,wind_num) = tempwt2.outputs(:,2);
    Pout2_all(:,wind_num) = tempwt2.outputs(:,3);
    Pitch2_all(:,wind_num) = tempwt2.states(:,1);
    Wr2_all(:,wind_num) = tempwt2.states(:,2);
    Wg2_all(:,wind_num) = tempwt2.states(:,3);
end
% Pref2_delay = zeros(10,100);
% Pref2_delaycount = zeros(10,100);
% Wind2_delay = zeros(10,100);
for wind_num=1:100
    index = find(abs(Pref2_all(1:end-1,wind_num)-Pref2_all(2:end,wind_num))==0);
    if ~isempty(index)
        timelength = length(index);
        Pref2_delay(1:timelength,wind_num) = index;
        flag=1;
        flag_num=1;
        for i=2:timelength
            if index(i)-index(i-1)==1
                flag_num = flag_num+1;
                
            else
                Pref2_delaycount(flag+1,wind_num)=flag_num;
                flag =flag+1;
                flag_num = 1;
            end
        end
        Pref2_delaycount(flag+1,wind_num)=flag_num;
        Pref2_delaycount(1,wind_num)=flag;
    end
    index = find(abs(Wind2_all(1:end-1,wind_num)-Wind2_all(2:end,wind_num))==0);
    if ~isempty(index)
         timelength = length(index);
        Wind2_delay(1:timelength,wind_num) = index;
        
        flag=1;
        flag_num=1;
        for i=2:timelength
            if index(i)-index(i-1)==1
                flag_num = flag_num+1;
                
            else
                Wind2_delaycount(flag+1,wind_num)=flag_num;
                flag =flag+1;
                flag_num = 1;
            end
        end
        Wind2_delaycount(flag+1,wind_num)=flag_num;
        Wind2_delaycount(1,wind_num)=flag;
        
    end
    index = find(abs(Tshaft2_all(1:end-1,wind_num)-Tshaft2_all(2:end,wind_num))==0);
    if ~isempty(index)
         timelength = length(index);
        Tshaft2_delay(1:timelength,wind_num) = index;
        
        flag=1;
        flag_num=1;
        for i=2:timelength
            if index(i)-index(i-1)==1
                flag_num = flag_num+1;
                
            else
                Tshaft2_delaycount(flag+1,wind_num)=flag_num;
                flag =flag+1;
                flag_num = 1;
            end
        end
        Tshaft2_delaycount(flag+1,wind_num)=flag_num;
        Tshaft2_delaycount(1,wind_num)=flag;
        
    end
    index = find(abs(Ft2_all(1:end-1,wind_num)-Ft2_all(2:end,wind_num))==0);
    if ~isempty(index)
         timelength = length(index);
        Ft2_delay(1:timelength,wind_num) = index;
        
        flag=1;
        flag_num=1;
        for i=2:timelength
            if index(i)-index(i-1)==1
                flag_num = flag_num+1;
                
            else
                Ft2_delaycount(flag+1,wind_num)=flag_num;
                flag =flag+1;
                flag_num = 1;
            end
        end
        Ft2_delaycount(flag+1,wind_num)=flag_num;
        Ft2_delaycount(1,wind_num)=flag;
        
    end
    
end

% delaycount = Wind2_delaycount;
delaycount = Pref2_delaycount;
% delaycount = Tshaft2_delaycount;
% delaycount = Ft2_delaycount;

mean(delaycount(1,:))
sum(delaycount(1,:))
% 统计每个数字的出现次数
data = delaycount(2:end,:);
counts = histcounts(data, 1:11); % 1到10的区间，使用11来确保10被包含

% 绘制直方图
figure;
bar(1:10, counts);
xlabel('数字');
ylabel('出现次数');
title('1到10各延迟出现次数直方图');
xlim([1, 10]);
grid on;

%% fit验证 原始
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

%
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
%% data 修复

%% fit验证 误差
datat = 1901;
X1 = Pref2(datat-1:end-1);
Y1 = Tshaft(datat:end);
% 假设 Y 是列向量, X1-X6 是特征矩阵中的列
Y2 = [Ft(datat:end)];  % 你的目标变量Y
X2 = [Wind2(datat:end) Pitch2(datat:end) Wr2(datat-1:end-1)];  % 特征变量矩阵

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

%
X1 = [Wind2(datat-2:end-2) Pref2(datat-2:end-2)];  % 特征变量矩阵
Y1 = Pitch2(datat:end);

X2 = [Wind2(datat-1:end-1) Pref2(datat-1:end-1)];
Y2 = Wr2(datat:end);

Y3 = [Ft(datat:end)];  % 你的目标变量Y
X3 = [Wind2(datat:end) Pitch2(datat:end) Wr2(datat-1:end-1)];  % 特征变量矩阵

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
%% fit验证 误差修正后
%修正数据
Fs = 1; 
Fc = 0.25;
% 设计低通滤波器
d = designfilt('lowpassfir', 'FilterOrder', 1, ...
                'CutoffFrequency', Fc, ...
                'SampleRate', Fs);

index = find(Pref2(1:end-1)-Pref2(2:end)==0);
for i=1:length(index)
    flag = index(i);
    if flag>10
        Pref2(flag) = mean(Pref2(flag-10:flag-1));
    else
        Pref2(flag) = mean(Pref2(1:flag-1));
    end
end
signal = Pref2;
% 使用滤波器
a=0.85;b=1-a;
Pref2 = b*filter(d, signal)+a*Pref2;

index = find(Wind2(1:end-1)-Wind2(2:end)==0);
for i=1:length(index)
    flag = index(i);
    if flag>10
        Wind2(flag) = mean(Wind2(flag-10:flag-1));
    else
        Wind2(flag) = mean(Wind2(1:flag-1));
    end
end
signal = Wind2;
% 使用滤波器
a=0.85;b=1-a;
Wind2 = b*filter(d, signal)+a*Wind2;
           
    

datat = 1901;
X1 = Pref2(datat-1:end-1);
Y1 = Tshaft(datat:end);
% 假设 Y 是列向量, X1-X6 是特征矩阵中的列
Y2 = [Ft(datat:end)];  % 你的目标变量Y
X2 = [Wind2(datat:end) Pitch2(datat:end) Wr2(datat-1:end-1)];  % 特征变量矩阵

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

%
X1 = [Wind2(datat-2:end-2) Pref2(datat-2:end-2)];  % 特征变量矩阵
Y1 = Pitch2(datat:end);

X2 = [Wind2(datat-1:end-1) Pref2(datat-1:end-1)];
Y2 = Wr2(datat:end);

Y3 = [Ft(datat:end)];  % 你的目标变量Y
X3 = [Wind2(datat:end) Pitch2(datat:end) Wr2(datat-1:end-1)];  % 特征变量矩阵

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
