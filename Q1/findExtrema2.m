% 定义寻找局部极大值和极小值的函数
function findExtrema2(current_index, turbine_index,raw_data)
    persistent n m Tmin_data Tmax_data minima_stack maxima_stack Ptrmax Ptrmin;
    persistent half_cycle_count full_cycle_count cycle_data;
    persistent last_count theta_b mm C N_total;
    global L_N D_F;
    if isempty(Tmin_data)
        n=100;m=100;
        Tmin_data = zeros(n, m);  % 存储极小值，二维矩阵
        Tmax_data = zeros(n, m);  % 存储极大值，二维矩阵
        minima_stack = zeros(n, m);  % 最小值堆栈，二维矩阵
        maxima_stack = zeros(n, m);  % 最大值堆栈，二维矩阵
        Ptrmax = zeros(1, m);  % 每台风机的最大值指针
        Ptrmin = zeros(1, m);  % 每台风机的最小值指针
        half_cycle_count = zeros(1, m);  % 每台风机的半循环计数器
        full_cycle_count = zeros(1, m);  % 每台风机的完整循环计数器
        last_count = zeros(1, m);
        cycle_data = zeros(n, 3, m);  % 三维矩阵，第三维表示风机的序号，存储半循环和全循环数据
        theta_b = 5e7;
        mm = 10;
        C = 9.77e70;
        N_total = 42565440.4361;
    end
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
%     rainflowCounting(current_index,turbine_index);

    % 常量定义
    
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
        L_p = (sum(S.^mm .* cycle_data(:, 1, turbine_index) ./ N_total))^(1/mm);
%         disp(['风机 ', num2str(turbine_index), ' 的 L_p = ', num2str(L_p)]);
        L_N(current_index, turbine_index) = L_p;
        
        % 计算N
        N = C ./ S.^mm;
        
        % 计算Df_p
        Df_p = sum(cycle_data(:, 1, turbine_index) ./ N);
        D_F(current_index, turbine_index) = Df_p;
        
%         disp(['风机 ', num2str(turbine_index), ' 的 Df_p = ', num2str(Df_p)]);
        last_count(turbine_index) = count;
    end
end

