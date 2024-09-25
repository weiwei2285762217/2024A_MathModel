
%%

%% Tshaft
% 定义全局变量，用于存储数据和计数
global  L_N D_F;
tic
% 初始化变量
n = 100;  % 数据长度为100
m = 100;  % 100台风机
L_N = zeros(n, m);  % 存储L_p值，二维矩阵
D_F = zeros(n, m);  % 存储Df_p值，二维矩阵
tempdata = Tshaft(1:n, 1:m); 

for tt = 1:100

    % 清除命令行窗口
%     clc;
    current_index = tt;
    % 确保 current_index 在有效范围内
    if current_index <= n && current_index > 0
        % 读取每台风机当前时刻的数据
        raw_data(current_index, :) = tempdata(current_index, :);
        
%         fprintf('读取第 %d 个数据\n', current_index);
        
        % 实时处理每台风机的数据
        for turbine_index = 1:m
            findExtrema2(current_index, turbine_index,raw_data);
        end
        
        % 更新索引
%         current_index = current_index + 1;
    else
        
    end
end

% % 定义定时器
% dataTimer = timer('ExecutionMode', 'fixedRate', 'Period', 0.001, 'TimerFcn', @readData);
% 
% % 启动定时器
% start(dataTimer);
toc
%%
figure()
plot(1:100,Tshaft(101,:))
hold on
plot(1:100,L_N(100,:))
legend('real','pre')

figure()
plot(1:100,Tshaft(102,:))
hold on
plot(1:100,D_F(100,:))
legend('real','pre')
ylim([0 2e-17])

figure()
plot(1:100,D_F(:,50))

%% Ft
% 定义全局变量，用于存储数据和计数
global  L_N D_F;
tic
% 初始化变量
n = 100;  % 数据长度为100
m = 100;  % 100台风机
L_N = zeros(n, m);  % 存储L_p值，二维矩阵
D_F = zeros(n, m);  % 存储Df_p值，二维矩阵
tempdata = Ft(1:n, 1:m); 

for tt = 1:100

    % 清除命令行窗口
%     clc;
    current_index = tt;
    % 确保 current_index 在有效范围内
    if current_index <= n && current_index > 0
        % 读取每台风机当前时刻的数据
        raw_data(current_index, :) = tempdata(current_index, :);
        
%         fprintf('读取第 %d 个数据\n', current_index);
        
        % 实时处理每台风机的数据
        for turbine_index = 1:m
           findExtrema2(current_index, turbine_index,raw_data);
        end
        
        % 更新索引
%         current_index = current_index + 1;
    else
        
    end
end

% % 定义定时器
% dataTimer = timer('ExecutionMode', 'fixedRate', 'Period', 0.001, 'TimerFcn', @readData);
% 
% % 启动定时器
% start(dataTimer);
toc
%%
figure()
yyaxis left
plot(1:100,Ft(101,:))
hold on
yyaxis right
plot(1:100,L_N(100,:))
legend('real','pre')

figure()
plot(1:100,Ft(102,:))
hold on
plot(1:100,D_F(100,:))
legend('real','pre')
%%
figure()
plot(1:100,D_F(:,50))

