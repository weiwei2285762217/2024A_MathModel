% 设置参数
mu = 0; % 均值
sigma = 0.25; % 标准差
n = 1000; % 生成随机数的数量

% 生成正态分布的随机数
data = normrnd(mu, sigma, [n, 1]);

% 将数据限制在-1到1之间
data(data < -1) = -1; % 将小于-1的值设为-1
data(data > 1) = 1;   % 将大于1的值设为1

% 绘制直方图
histogram(data, 'Normalization', 'probability');
xlabel('Value');
ylabel('Probability');
title('Histogram of Normal Distribution Random Numbers between -1 and 1');
