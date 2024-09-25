% 1. 读取Excel数据
temperature_data = [-2 1 -3 5 -1 3 -4 4 -3 1 -2 3 2 6];
temperature_data = Tshaft(1:100,100);
L = Tshaft(101,1);
raw_data = temperature_data;  % 假设表格中存储的是温度数据
n = length(raw_data);  % 数据的长度

% 初始化变量
Tmin_data = [];  % 用于存储极小值
Tmax_data = [];  % 用于存储极大值

global minima_stack maxima_stack Ptrmax Ptrmin half_cycle_count full_cycle_count;
% 2. 进行极值处理：找到数据中的局部极大值和极小值
for i = 2:n-1
    if raw_data(i) > raw_data(i-1) && raw_data(i) > raw_data(i+1)
        % 如果当前值大于前后值，则为极大值
        Tmax_data = [Tmax_data; raw_data(i)];
    elseif raw_data(i) < raw_data(i-1) && raw_data(i) < raw_data(i+1)
        % 如果当前值小于前后值，则为极小值
        Tmin_data = [Tmin_data; raw_data(i)];
    end
end
%% 

% 检查第一个和最后一个数据点
if raw_data(1) > raw_data(2)
    Tmax_data = [Tmax_data; raw_data(1)];  % 如果第一个数据是极大值
elseif raw_data(1) < raw_data(2)
    Tmin_data = [Tmin_data; raw_data(1)];  % 如果第一个数据是极小值
end

if raw_data(n) > raw_data(n-1)
    Tmax_data = [Tmax_data; raw_data(n)];  % 如果最后一个数据是极大值
elseif raw_data(n) < raw_data(n-1)
    Tmin_data = [Tmin_data; raw_data(n)];  % 如果最后一个数据是极小值
end
%% 

% 初始化堆栈和指针
minima_stack = [];  % 存储最小值的堆栈
maxima_stack = [];  % 存储最大值的堆栈
Ptrmax = 0;  % 最大值指针
Ptrmin = 0;  % 最小值指针

% 初始化循环计数器
half_cycle_count = 0;    % 半循环计数器
full_cycle_count = 0;    % 完整循环计数器



% 3. 将极大值和极小值数据传入雨流计数流程
Tmin_data_length = length(Tmin_data);
Tmax_data_length = length(Tmax_data);
min_length = min(Tmin_data_length, Tmax_data_length);

for i = 1:min_length
    Tmin_new = Tmin_data(i);  % 从极小值数据中读取当前的Tmin
    Tmax_new = Tmax_data(i);  % 从极大值数据中读取当前的Tmax
    
    % 处理最小值
    if isempty(minima_stack) || Tmin_new < minima_stack(end)
        % 如果新的Tmin小于当前堆栈中的最小值，或堆栈为空
        minima_stack = [minima_stack, Tmin_new];  % 新的最小值入栈
        Ptrmin = Ptrmin + 1;  % 更新最小值指针
        
        % 处理最大值
        if isempty(maxima_stack)
            continue;  % 如果没有最大值，等待最大值更新
        else
            % 检查是否存在完整或半循环
            recursiveRainflow();  % 调用递归方法处理堆栈中的值
        end
    elseif Tmin_new >= minima_stack(end)
        % 如果新的Tmin比堆栈中的最小值大，则将新值保存到最小值堆栈
        minima_stack = [minima_stack, Tmin_new];
        Ptrmin = Ptrmin + 1;  % 更新最小值指针
    end
    
    % 处理最大值
    
    if isempty(maxima_stack) || Tmax_new > maxima_stack(end)
        maxima_stack = [maxima_stack, Tmax_new];  % 新的最大值入栈
        Ptrmax = Ptrmax + 1;  % 更新最大值指针
        Ptrmax
        maxima_stack
    end
end

% 输出半循环和完整循环的次数
disp(['半循环次数: ', num2str(half_cycle_count)]);
disp(['完整循环次数: ', num2str(full_cycle_count)]);


% 定义递归函数处理最小值和最大值
function recursiveRainflow()
    global minima_stack maxima_stack Ptrmax Ptrmin half_cycle_count full_cycle_count;
    
    while ~isempty(minima_stack) && ~isempty(maxima_stack)
        % 取得栈顶最小值和最大值
        Tmin_old = minima_stack(end);
        Tmax_old = maxima_stack(end);
        
        % 如果指针Ptrmax为1，表示存在一个最大值 -> 半循环
        if Ptrmax == 1
            delta_T = abs(Tmax_old - Tmin_old);  % 计算半循环的温度变化
            half_cycle_count = half_cycle_count + 1;  % 增加半循环计数
            
            % 移除最小值
            minima_stack(end) = [];  % 移除旧的最小值
            Ptrmin = Ptrmin - 1;  % 更新最小值指针
        elseif Ptrmax > 1
            % 如果Ptrmax>1，表示存在多个最大值 -> 完整循环
            Tmax_new = maxima_stack(end-1);  % 取得第二个最大值
            delta_T = abs(Tmax_new - Tmin_old);  % 计算完整循环的温度变化
            full_cycle_count = full_cycle_count + 1;  % 增加完整循环计数
            
            % 移除最小值和新的最大值
            minima_stack(end) = [];  % 移除旧的最小值
            maxima_stack(end) = [];  % 移除新的最大值
            Ptrmin = Ptrmin - 1;  % 更新最小值指针
            Ptrmax = Ptrmax - 1;  % 更新最大值指针
        end
        
        % 如果 Ptrmin 指针大于等于1，递归地重复整个操作
        if Ptrmin >= 1
            recursiveRainflow();  % 递归调用
        end
    end
end
