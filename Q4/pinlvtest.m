%% 频率相干性
% 假设 model1 和 model2 是你的两个信号数据
% model1 = ... (第一个型号的信号数据)
% model2 = ... (第二个型号的信号数据)

% 采样频率
Fs = 1; % 假设采样频率为 1 Hz
model1 = Wind-mean(Wind);
model2 = Wind2-mean(Wind2);
% 计算信号的长度
N = length(model1);

% 使用 pwelch 计算功率谱密度和交叉谱密度
[Px, ~] = pwelch(model1, [], [], [], Fs); % 第一个信号的功率谱密度
[Py, ~] = pwelch(model2, [], [], [], Fs); % 第二个信号的功率谱密度
[Pxy, F] = cpsd(model1, model2, [], [], [], Fs); % 交叉谱密度

% 计算相干函数
Cxy = abs(Pxy).^2 ./ (Px .* Py); % 相干函数

% 绘制相干性图
figure;
plot(F, Cxy, 'LineWidth', 1.5);
title('频率相干性');
xlabel('频率 (Hz)');
ylabel('相干性');
xlim([0 Fs/2]); % 限制频率范围
grid on;

% 显示最大相干性
disp(['最大相干性: ', num2str(max(Cxy))]);


%% 原始
% 采样频率
Fs = 1; % 这里假设采样频率为 1 Hz
signal = Wind-mean(Wind);
% 计算信号的长度
N = length(signal);

% 进行快速傅里叶变换
Y = fft(signal);

% 计算频率轴
f = (0:N-1)*(Fs/N); % 频率向量

% 计算双侧幅度谱
P2 = abs(Y/N); % 双侧幅度谱

% 计算单侧幅度谱
P1 = P2(1:N/2+1); % 单侧幅度谱
P1(2:end-1) = 2*P1(2:end-1); % 去掉重复的部分

% 绘图
figure;
plot(f(1:N/2+1), P1);
title('单向频率图');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;

%% 扰动
Fs = 1; % 这里假设采样频率为 1 Hz
signal = Wind2-mean(Wind2);
% 计算信号的长度
N = length(signal);

% 进行快速傅里叶变换
Y = fft(signal);

% 计算频率轴
f = (0:N-1)*(Fs/N); % 频率向量

% 计算双侧幅度谱
P2 = abs(Y/N); % 双侧幅度谱

% 计算单侧幅度谱
P1 = P2(1:N/2+1); % 单侧幅度谱
P1(2:end-1) = 2*P1(2:end-1); % 去掉重复的部分

% 绘图
figure;
plot(f(1:N/2+1), P1);
title('单向频率图');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;
