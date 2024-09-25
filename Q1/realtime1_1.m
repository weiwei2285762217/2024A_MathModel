
%%
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
%%
figure()
plot(1:100,D_F(:,50))
%%
% 定义全局变量，用于存储数据和计数
global raw_data n current_index temperature_data Tmin_data Tmax_data minima_stack maxima_stack Ptrmax Ptrmin half_cycle_count full_cycle_count cycle_data L_N D_F m;
global last_count
tic
% 初始化变量
n = 100;  % 数据长度为100
m = 100;  % 100台风机
temperature_data = Tshaft(1:n, 1:m);  % 假设Tshaft是n*m的温度数据矩阵
raw_data = zeros(n, m);  % 存储读取的温度数据
Tmin_data = zeros(n, m);  % 存储极小值，二维矩阵
Tmax_data = zeros(n, m);  % 存储极大值，二维矩阵
minima_stack = zeros(n, m);  % 最小值堆栈，二维矩阵
maxima_stack = zeros(n, m);  % 最大值堆栈，二维矩阵
Ptrmax = zeros(1, m);  % 每台风机的最大值指针
Ptrmin = zeros(1, m);  % 每台风机的最小值指针
half_cycle_count = zeros(1, m);  % 每台风机的半循环计数器
full_cycle_count = zeros(1, m);  % 每台风机的完整循环计数器
last_count = zeros(1, m);
current_index = 1;  % 当前读取数据的索引
cycle_data = zeros(n, 3, m);  % 三维矩阵，第三维表示风机的序号，存储半循环和全循环数据
L_N = zeros(n, m);  % 存储L_p值，二维矩阵
D_F = zeros(n, m);  % 存储Df_p值，二维矩阵


for tt = 1:100
    global raw_data current_index temperature_data n m;

    % 清除命令行窗口
%     clc;
    current_index = tt;
    % 确保 current_index 在有效范围内
    if current_index <= n && current_index > 0
        % 读取每台风机当前时刻的数据
        raw_data(current_index, :) = temperature_data(current_index, :);
        
%         fprintf('读取第 %d 个数据\n', current_index);
        
        % 实时处理每台风机的数据
        for turbine_index = 1:m
            findExtrema(current_index, turbine_index);
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
% 定义读取数据的回调函数
function readData(~, ~)
    global raw_data current_index temperature_data n m;

    % 清除命令行窗口
%     clc;

    % 确保 current_index 在有效范围内
    if current_index <= n && current_index > 0
        % 读取每台风机当前时刻的数据
        raw_data(current_index, :) = temperature_data(current_index, :);
        
        fprintf('读取第 %d 个数据\n', current_index);
        
        % 实时处理每台风机的数据
        for turbine_index = 1:m
            findExtrema(current_index, turbine_index);
        end
        
        % 更新索引
        current_index = current_index + 1;
    else
        % 当所有数据读取完毕时，停止定时器
        disp('所有数据已读取完毕。');
    end
end

% 定义寻找局部极大值和极小值的函数
function findExtrema(current_index, turbine_index)
    global raw_data Tmin_data Tmax_data minima_stack maxima_stack Ptrmax Ptrmin;

    % 检查 current_index 和 turbine_index 是否为有效索引
    if current_index < 3 || current_index > size(raw_data, 1) || turbine_index > size(raw_data, 2)
        
        return;
    end
    
    
    % 检查第 current_index-1 个数据是否为局部极大值或极小值
    i = current_index - 1;  % 当前检查的索引为当前已读取数据的倒数第二个
    
    if raw_data(i, turbine_index) > raw_data(i-1, turbine_index) && raw_data(i, turbine_index) > raw_data(i+1, turbine_index)
        % 局部极大值
        Ptrmax(turbine_index) = Ptrmax(turbine_index) + 1;
        Tmax_data(Ptrmax(turbine_index), turbine_index) = raw_data(i, turbine_index);
        maxima_stack(Ptrmax(turbine_index), turbine_index) = raw_data(i, turbine_index);
    elseif raw_data(i, turbine_index) < raw_data(i-1, turbine_index) && raw_data(i, turbine_index) < raw_data(i+1, turbine_index)
        % 局部极小值
        Ptrmin(turbine_index) = Ptrmin(turbine_index) + 1;
        Tmin_data(Ptrmin(turbine_index), turbine_index) = raw_data(i, turbine_index);
        minima_stack(Ptrmin(turbine_index), turbine_index) = raw_data(i, turbine_index);
    end
    
    % 进行雨流计数并计算
    rainflowCounting(current_index,turbine_index);
end

% 定义雨流计数函数并在每次循环后进行计算
function rainflowCounting(current_index,turbine_index)
    global n minima_stack maxima_stack Ptrmax Ptrmin half_cycle_count full_cycle_count cycle_data L_N D_F;
    global last_count;
    % 常量定义
    theta_b = 5e7;
    m = 10;
    C = 9.77e70;
    N_total = 42565440.4361;
    
    % 确保 Ptrmax 和 Ptrmin 是有效索引
    if Ptrmin(turbine_index) < 1 || Ptrmax(turbine_index) < 1
        L_N(current_index, turbine_index) = L_N(current_index-1, turbine_index);
        D_F(current_index, turbine_index) = D_F(current_index-1, turbine_index);
        return;
    end
    
    while Ptrmin(turbine_index) >= 1 && Ptrmax(turbine_index) >= 1
        Tmin_old = minima_stack(Ptrmin(turbine_index), turbine_index);
        Tmax_old = maxima_stack(Ptrmax(turbine_index), turbine_index);
        
        % 如果只有一个最大值，半个周期
        if Ptrmax(turbine_index) == 1
            delta_T = abs(Tmax_old - Tmin_old);
            bar_T = (Tmax_old + Tmin_old) / 2;
            half_cycle_count(turbine_index) = half_cycle_count(turbine_index) + 1;
            
            % 将半循环数据存储到 cycle_data 矩阵，第三维为风机序号
            cycle_data(half_cycle_count(turbine_index)+full_cycle_count(turbine_index), :, turbine_index) = [0.5, delta_T, bar_T];
            
            % 移除最小值
            minima_stack(Ptrmin(turbine_index), turbine_index) = 0;
            Ptrmin(turbine_index) = Ptrmin(turbine_index) - 1;
        elseif Ptrmax(turbine_index) > 1
            % 完整周期
            Tmax_new = maxima_stack(Ptrmax(turbine_index)-1, turbine_index);
            delta_T = abs(Tmax_new - Tmin_old);
            bar_T = (Tmax_new + Tmin_old) / 2;
            full_cycle_count(turbine_index) = full_cycle_count(turbine_index) + 1;
            
            % 将完整循环数据存储到 cycle_data 矩阵
            cycle_data(full_cycle_count(turbine_index)+half_cycle_count(turbine_index), :, turbine_index) = [1, delta_T, bar_T];
            
            % 移除最小值和新最大值
            minima_stack(Ptrmin(turbine_index), turbine_index) = 0;
            maxima_stack(Ptrmax(turbine_index), turbine_index) = 0;
            Ptrmin(turbine_index) = Ptrmin(turbine_index) - 1;
            Ptrmax(turbine_index) = Ptrmax(turbine_index) - 1;
        end
        
    end
    
    % 输出半循环和完整循环的次数
%     disp(['风机 ', num2str(turbine_index), ' 的半循环次数: ', num2str(half_cycle_count(turbine_index))]);
%     disp(['风机 ', num2str(turbine_index), ' 的完整循环次数: ', num2str(full_cycle_count(turbine_index))]);

    % ---- 执行计算 ----
    count = half_cycle_count(turbine_index)+full_cycle_count(turbine_index);
    if last_count(turbine_index)<count
%     if ~isempty(cycle_data)
        % 计算S
        S = cycle_data(:, 2, turbine_index) ./ (1 - cycle_data(:, 3, turbine_index) ./ theta_b);
        
        % 计算L_p
        L_p = (sum(S.^m .* cycle_data(:, 1, turbine_index) ./ N_total))^(1/m);
%         disp(['风机 ', num2str(turbine_index), ' 的 L_p = ', num2str(L_p)]);
        L_N(current_index, turbine_index) = L_p;
        
        % 计算N
        N = C ./ S.^m;
        
        % 计算Df_p
        Df_p = sum(cycle_data(:, 1, turbine_index) ./ N);
        D_F(current_index, turbine_index) = Df_p;
        
%         disp(['风机 ', num2str(turbine_index), ' 的 Df_p = ', num2str(Df_p)]);
        last_count(turbine_index) = count;
    end
end

