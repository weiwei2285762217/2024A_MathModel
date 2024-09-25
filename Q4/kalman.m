% 生成模拟数据
N = 100; % 数据点数
x = (1:N)'; % 自变量 x
a_true = 2; % 真实的斜率 a
b_true = 5; % 真实的截距 b
y = a_true * x + b_true + randn(N, 1); % 加入噪声后的 y

% 卡尔曼滤波参数
A = eye(2);      % 状态转移矩阵，假设 a 和 b 不随时间变化
Q = 1e-5 * eye(2); % 过程噪声协方差矩阵，假设小噪声
R = 1;           % 观测噪声方差

% 初始化状态
theta_est = [0; 0]; % 初始估计的 a 和 b
P = eye(2);        % 初始协方差矩阵

% 存储估计的 a 和 b
theta_estimates = zeros(N, 2);

% 卡尔曼滤波过程
for k = 1:N
    % 当前观测矩阵 H_k
    H_k = [x(k), 1]; 
    
    % 预测步骤
    theta_pred = A * theta_est;
    P_pred = A * P * A' + Q;
    
    % 计算卡尔曼增益
    K = P_pred * H_k' / (H_k * P_pred * H_k' + R);
    
    % 更新状态估计
    theta_est = theta_pred + K * (y(k) - H_k * theta_pred);
    
    % 更新协方差矩阵
    P = (eye(2) - K * H_k) * P_pred;
    
    % 存储当前估计的 a 和 b
    theta_estimates(k, :) = theta_est';
end

% 绘制估计的 a 和 b
figure;
subplot(2,1,1);
plot(1:N, theta_estimates(:, 1), 'b-', 'LineWidth', 1.5);
hold on;
plot(1:N, a_true * ones(N, 1), 'r--', 'LineWidth', 1.5);
legend('估计的 a', '真实的 a');
title('卡尔曼滤波对 a 的估计');

subplot(2,1,2);
plot(1:N, theta_estimates(:, 2), 'b-', 'LineWidth', 1.5);
hold on;
plot(1:N, b_true * ones(N, 1), 'r--', 'LineWidth', 1.5);
legend('估计的 b', '真实的 b');
title('卡尔曼滤波对 b 的估计');
