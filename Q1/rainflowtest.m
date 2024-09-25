%% three point
allL = zeros(1,100);
allDf = zeros(1,100);
tic
for i=1:100
    Ft1 = Tshaft(1:100,i);
    load_time_series = Ft1;  % 原始载荷数据
    
    c = three_point_rainflow(load_time_series);
%     c = three_point_rainflow(load_time_series);
    %goodman
    theta_b = 5e7;
    S = c(:,2)./(1-c(:,3)./theta_b);
    %
    m = 10;
    N = 42565440.4361;
%     L
    allL(1,i) = (sum(S.^m.*c(:,1)./N))^(1/m);
    %
    C = 9.77e70;
    N = (C./S.^m);
%     Df
    allDf(1,i) = sum(c(:,1)./N);
end
L = Ft(101,:);
Df = Ft(102,:);
toc
figure()

plot(1:100,allL)
hold on
plot(1:100,L)
title('L')
legend('pre','real')
figure()

plot(1:100,allDf)
hold on
plot(1:100,Df)
legend('pre','real')
title('Df')
%% rainflow--all
allL = zeros(1,100);
allDf = zeros(1,100);
for i=1:100
    Ft1 = Ft(1:100,i);
    load_time_series = Ft1;  % 原始载荷数据
    time = 1:length(load_time_series);  % 时间序列

    % 提取波峰
    [peaks, peak_locs] = findpeaks(load_time_series, time);

    % 提取波谷 (通过对数据取反找波峰)
    [valleys, valley_locs] = findpeaks(-load_time_series, time);
    valleys = -valleys;  % 将波谷恢复为原始值

    % 合并波峰和波谷
    extrema_values = [peaks; valleys];
    extrema_locs = [peak_locs'; valley_locs'];

    % 按时间排序波峰波谷
    [extrema_locs_sorted, sort_idx] = sort(extrema_locs);
    extrema_values_sorted = extrema_values(sort_idx);
    
    c = rainflow(extrema_values_sorted,extrema_locs_sorted);

    %goodman
    theta_b = 5e7;
    S = c(:,2)./(1-c(:,3)./theta_b);
    %
    m = 10;
    N = 42565440.4361;
%     L
    allL(1,i) = (sum(S.^m.*c(:,1)./N))^(1/m);
    %
    C = 9.77e70;
    N = (C./S.^m);
%     Df
    allDf(1,i) = sum(c(:,1)./N);
end
L = Ft(101,:);
Df = Ft(102,:);
figure()

plot(1:100,allL)
hold on
plot(1:100,L)
title('L')
legend('pre','real')
figure()

plot(1:100,allDf)
hold on
plot(1:100,Df)
legend('pre','real')
title('Df')
%% rainflow--all
allL = zeros(1,100);
allDf = zeros(1,100);
for i=1:100
    Ft1 = Tshaft(1:100,i);
    load_time_series = Ft1;  % 原始载荷数据
    time = 1:length(load_time_series);  % 时间序列

    % 提取波峰
    [peaks, peak_locs] = findpeaks(load_time_series, time);

    % 提取波谷 (通过对数据取反找波峰)
    [valleys, valley_locs] = findpeaks(-load_time_series, time);
    valleys = -valleys;  % 将波谷恢复为原始值

    % 合并波峰和波谷
    extrema_values = [peaks; valleys];
    extrema_locs = [peak_locs'; valley_locs'];

    % 按时间排序波峰波谷
    [extrema_locs_sorted, sort_idx] = sort(extrema_locs);
    extrema_values_sorted = extrema_values(sort_idx);
    
    c = rainflow(extrema_values_sorted,extrema_locs_sorted);

    %goodman
    theta_b = 5e7;
    S = c(:,2)./(1-c(:,3)./theta_b);
    %
    m = 10;
    N = 42565440.4361;
%     L
    allL(1,i) = (sum(S.^m.*c(:,1)./N))^(1/m);
    %
    C = 9.77e70;
    N = (C./S.^m);
%     Df
    allDf(1,i) = sum(c(:,1)./N);
end
L = Tshaft(101,:);
Df = Tshaft(102,:);
figure()

plot(1:100,allL)
hold on
plot(1:100,L)
title('L')
legend('pre','real')
figure()

plot(1:100,allDf)
hold on
plot(1:100,Df)
legend('pre','real')
title('Df')