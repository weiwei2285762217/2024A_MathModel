
%% 采样频率
Fs = 1; % 1 Hz

% 截止频率
Fc = 0.4; % 0.25 Hz

% 设计低通滤波器
d = designfilt('lowpassfir', 'FilterOrder', 1, ...
                'CutoffFrequency', Fc, ...
                'SampleRate', Fs);
% 示例信号
wt_num = 100;
index = find(Wind(1:end-1,wt_num)-Wind(2:end,wt_num)==0);
for i=1:length(index)
    flag = index(i);
    Wind(flag,wt_num) = mean(Wind(flag-10:flag-1,wt_num));
end

signal = Wind(:,wt_num);
t = 1:length(Wind(:,wt_num));
% 使用滤波器
filtered_signal = filter(d, signal);

% 绘图
figure;

a=0.5;b=1-a;
plot(t, Wind(:,wt_num));
hold on
plot(t, b*filtered_signal+a*Wind(:,wt_num));
plot(t,Wind3(:,wt_num))
xlabel('时间 (s)');
ylabel('幅值');
legend('噪声','滤波','真实')
xlim([0 100])

real = Wind3(:,wt_num);pre = Wind(:,wt_num);
mse = sum((real - pre).^2) / length(real)
mape = sum(abs((real - pre) ./ real)) / length(real) * 100



real = Wind3(:,wt_num);pre = b*filtered_signal+a*Wind(:,wt_num);
mse = sum((real - pre).^2) / length(real)
mape = sum(abs((real - pre) ./ real)) / length(real) * 100
%% 采样频率
Fs = 1; % 1 Hz

% 截止频率
Fc = 0.25; % 0.25 Hz

% 设计低通滤波器
d = designfilt('lowpassfir', 'FilterOrder', 1, ...
                'CutoffFrequency', Fc, ...
                'SampleRate', Fs);
% 示例信号
wt_num = 1;
index = find(Pref(1:end-1,wt_num)-Pref(2:end,wt_num)==0);
for i=1:length(index)
    flag = index(i);
    Pref(flag,wt_num) = mean(Pref(flag-10:flag-1,wt_num));
end

signal = Pref(:,wt_num);
t = 1:length(Pref(:,wt_num));
% 使用滤波器
filtered_signal = filter(d, signal);

% 绘图
figure;

a=0.85;b=1-a;
plot(t, Pref(:,wt_num));
hold on
plot(t, b*filtered_signal+a*Pref(:,wt_num));
plot(t,Pref3(:,wt_num))
xlabel('时间 (s)');
ylabel('幅值');
legend('噪声','滤波','真实')
xlim([0 100])

real = Pref3(:,wt_num);pre = Pref(:,wt_num);
mse = sum((real - pre).^2) / length(real)
mape = sum(abs((real - pre) ./ real)) / length(real) * 100

real = Pref3(:,wt_num);pre = b*filtered_signal+a*Pref(:,wt_num);
mse = sum((real - pre).^2) / length(real)
mape = sum(abs((real - pre) ./ real)) / length(real) * 100