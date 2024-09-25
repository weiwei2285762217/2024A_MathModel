function c = three_point_rainflow(stress)
    % three_point_rainflow 使用三点式雨流计数法计算应力循环
    % 输入:
    %   stress - 应力时间序列 (一维数组)
    % 输出:
    %   count - 每个循环的计数
    %   ranges - 每个循环的幅值范围

    % 初始化
%     stress = extrema_values_sorted;  % 确保为列向量
    count = [];
    ranges = [];
    means = [];
    countnum=0;
    count1 = [];
    count2 = [];
    
    time = 1:length(stress);  % 时间序列
    
    % 提取波峰
    [peaks, peak_locs] = findpeaks(stress, time);

    % 提取波谷 (通过对数据取反找波峰)
    [valleys, valley_locs] = findpeaks(-stress, time);
    valleys = -valleys;  % 将波谷恢复为原始值

    % 合并波峰和波谷
    extrema_values = [peaks; valleys];
    extrema_locs = [peak_locs'; valley_locs'];

    % 按时间排序波峰波谷
    [extrema_locs_sorted, sort_idx] = sort(extrema_locs);
    extrema_values_sorted = extrema_values(sort_idx);
    stress_old = stress;
    stress = extrema_values_sorted;
    while length(stress) >= 3
        found_cycle = false;  % 标记是否找到循环

        % 依次考察三个连续的点
        for i = 1:length(stress) - 2
            S1 = stress(i);
            S2 = stress(i + 1);
            S3 = stress(i + 2);
            
            % 计算应力变化量
            deltaS1 = abs(S1 - S2);
            deltaS2 = abs(S2 - S3);
            
            % 判断是否形成有效循环
            if deltaS1 <= deltaS2
                % 形成循环，计算幅值并计数
                mean_value = (S1 + S2) / 2;  % 计算均值
                range = (max(S1, S2) - min(S1, S2));
                count(end + 1) = 1; % 循环次数
                ranges(end + 1) = range;
                means(end + 1) = mean_value;  % 存储均值
                 
                % 移除已识别的循环点
                count1(end+1) = find(stress_old==stress(i));
                count2(end+1) = find(stress_old==stress(i+1));
                stress(i:i + 1) = [];  % 移除 S1, S2
                
                found_cycle = true;  % 标记找到循环
                countnum = countnum+2;
                break;  % 重新开始
            end
        end
        
        % 如果没有找到循环，终止循环
        if ~found_cycle
            break;
        end
    end
    c= [count' ranges' means' count1' count2'];
    % 显示结果
%     disp('循环计数:');
%     disp(count);
%     disp('循环幅值:');
%     disp(ranges);
end


