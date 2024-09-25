% 定义全局变量，用于存储数据和计数
global raw_data n current_index temperature_data Tmin_data Tmax_data minima_stack maxima_stack Ptrmax Ptrmin half_cycle_count full_cycle_count;

% 初始化变量
n = 100;  % 数据长度为100
temperature_data = Tshaft(1:n, 100);  % 假设Tshaft的第100列是温度数据
raw_data = zeros(1, n);  % 存储读取的数据，初始化为0
Tmin_data = [];  % 存储极小值
Tmax_data = [];  % 存储极大值
minima_stack = [];  % 最小值堆栈
maxima_stack = [];  % 最大值堆栈
Ptrmax = 0;  % 最大值指针
Ptrmin = 0;  % 最小值指针
half_cycle_count = 0;    % 半循环计数器
full_cycle_count = 0;    % 完整循环计数器
current_index = 1;  % 当前读取数据的索引

% 定义定时器
dataTimer = timer('ExecutionMode', 'fixedRate', 'Period', 0.01, 'TimerFcn', @readData);

% 启动定时器
start(dataTimer);

% 定义读取数据的回调函数
function readData(~, ~)
    global raw_data current_index temperature_data n;

    % 清除命令行窗口
    clc;

    if current_index <= n
        % 每秒读取一个数据
        raw_data(current_index) = temperature_data(current_index);
        fprintf('读取第 %d 个数据: %f\n', current_index, raw_data(current_index));
        
        % 实时处理当前已读取的数据，寻找极大值和极小值
        findExtrema(current_index);
        
        % 更新索引
        current_index = current_index + 1;
    else
        % 当所有数据读取完毕时，停止定时器
        stop(timerfind);  % 停止所有定时器
        disp('所有数据已读取完毕。');
    end
end

% 定义寻找局部极大值和极小值的函数
function findExtrema(current_index)
    global raw_data Tmin_data Tmax_data minima_stack maxima_stack Ptrmax Ptrmin;

    % 如果已读取数据少于3个，则无法计算极值，直接返回
    if current_index < 3
        return;
    end
    
    % 检查第 current_index-1 个数据是否为局部极大值或极小值
    i = current_index - 1;  % 当前检查的索引为当前已读取数据的倒数第二个
    
    if raw_data(i) > raw_data(i-1) && raw_data(i) > raw_data(i+1)
        % 局部极大值
        Tmax_data = [Tmax_data; raw_data(i)];
        maxima_stack = [maxima_stack, raw_data(i)];
        Ptrmax = Ptrmax + 1;
    elseif raw_data(i) < raw_data(i-1) && raw_data(i) < raw_data(i+1)
        % 局部极小值
        Tmin_data = [Tmin_data; raw_data(i)];
        minima_stack = [minima_stack, raw_data(i)];
        Ptrmin = Ptrmin + 1;
    end
    
    % 进行雨流计数
    rainflowCounting();
    
    % 输出当前的极大值和极小值
    disp('当前局部极小值:');
    disp(Tmin_data);
    disp('当前局部极大值:');
    disp(Tmax_data);
end

% 定义雨流计数函数
function rainflowCounting()
    global minima_stack maxima_stack Ptrmax Ptrmin half_cycle_count full_cycle_count;
    
    while ~isempty(minima_stack) && ~isempty(maxima_stack)
        Tmin_old = minima_stack(end);
        Tmax_old = maxima_stack(end);
        
        % 如果只有一个最大值，半个周期
        if Ptrmax == 1
            delta_T = abs(Tmax_old - Tmin_old);
            half_cycle_count = half_cycle_count + 1;
            disp(['半循环: ΔT = ', num2str(delta_T)]);
            
            % 移除最小值
            minima_stack(end) = [];
            Ptrmin = Ptrmin - 1;
        elseif Ptrmax > 1
            % 完整周期
            Tmax_new = maxima_stack(end-1);
            delta_T = abs(Tmax_new - Tmin_old);
            full_cycle_count = full_cycle_count + 1;
            disp(['完整循环: ΔT = ', num2str(delta_T)]);
            
            % 移除最小值和新最大值
            minima_stack(end) = [];
            maxima_stack(end) = [];
            Ptrmin = Ptrmin - 1;
            Ptrmax = Ptrmax - 1;
        end
        
        % 检查是否需要递归处理
        if Ptrmin >= 1
            rainflowCounting();
        end
    end
    
    % 输出半循环和完整循环的次数
    disp(['半循环次数: ', num2str(half_cycle_count)]);
    disp(['完整循环次数: ', num2str(full_cycle_count)]);
end
