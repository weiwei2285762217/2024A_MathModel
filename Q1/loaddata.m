%% Tq
clear;clc;
opts = spreadsheetImportOptions("NumVariables", 100);

% 指定工作表和范围
opts.Sheet = "Sheet1";
opts.DataRange = "B2:CW103";

% 指定列名称和类型
opts.VariableNames = ["WT1", "WT2", "WT3", "WT4", "WT5", "WT6", "WT7", "WT8", "WT9", "WT10", "WT11", "WT12", "WT13", "WT14", "WT15", "WT16", "WT17", "WT18", "WT19", "WT20", "WT21", "WT22", "WT23", "WT24", "WT25", "WT26", "WT27", "WT28", "WT29", "WT30", "WT31", "WT32", "WT33", "WT34", "WT35", "WT36", "WT37", "WT38", "WT39", "WT40", "WT41", "WT42", "WT43", "WT44", "WT45", "WT46", "WT47", "WT48", "WT49", "WT50", "WT51", "WT52", "WT53", "WT54", "WT55", "WT56", "WT57", "WT58", "WT59", "WT60", "WT61", "WT62", "WT63", "WT64", "WT65", "WT66", "WT67", "WT68", "WT69", "WT70", "WT71", "WT72", "WT73", "WT74", "WT75", "WT76", "WT77", "WT78", "WT79", "WT80", "WT81", "WT82", "WT83", "WT84", "WT85", "WT86", "WT87", "WT88", "WT89", "WT90", "WT91", "WT92", "WT93", "WT94", "WT95", "WT96", "WT97", "WT98", "WT99", "WT100"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% 导入数据
Tshaft = readtable("H:\DATA\A\Q1\Tq.xlsx", opts, "UseExcel", false);

Tshaft = table2array(Tshaft);

clear opts
%Ft
opts = spreadsheetImportOptions("NumVariables", 100);

% 指定工作表和范围
opts.Sheet = "Sheet1";
opts.DataRange = "B2:CW103";

% 指定列名称和类型
opts.VariableNames = ["WT1", "WT2", "WT3", "WT4", "WT5", "WT6", "WT7", "WT8", "WT9", "WT10", "WT11", "WT12", "WT13", "WT14", "WT15", "WT16", "WT17", "WT18", "WT19", "WT20", "WT21", "WT22", "WT23", "WT24", "WT25", "WT26", "WT27", "WT28", "WT29", "WT30", "WT31", "WT32", "WT33", "WT34", "WT35", "WT36", "WT37", "WT38", "WT39", "WT40", "WT41", "WT42", "WT43", "WT44", "WT45", "WT46", "WT47", "WT48", "WT49", "WT50", "WT51", "WT52", "WT53", "WT54", "WT55", "WT56", "WT57", "WT58", "WT59", "WT60", "WT61", "WT62", "WT63", "WT64", "WT65", "WT66", "WT67", "WT68", "WT69", "WT70", "WT71", "WT72", "WT73", "WT74", "WT75", "WT76", "WT77", "WT78", "WT79", "WT80", "WT81", "WT82", "WT83", "WT84", "WT85", "WT86", "WT87", "WT88", "WT89", "WT90", "WT91", "WT92", "WT93", "WT94", "WT95", "WT96", "WT97", "WT98", "WT99", "WT100"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% 导入数据
Ft = readtable("H:\DATA\A\Q1\Ft.xlsx", opts, "UseExcel", false);
Ft = table2array(Ft);
clear opts

%% dataview
n=50;
time = 1:1:100;
Ft1 = Ft(1:100,n);
L = Ft(101,n);
Df = Ft(102,n);

Tshaft1 = Tshaft(1:100,n);
% L = Ft(101,n);
% Df = Ft(102,n);
%% 频域
% Fs = 1;
% X = Ft1-mean(Ft1);
% % X = highpass(Ft1, 0.1, Fs);  % fs 是采样频率
% len = length(X);
% Y = fft(X);
% P2 = abs(Y/len);
% P1 = P2(1:len/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% figure()
% plot(Fs/len*(0:len/2),P1,"LineWidth",3)
% title("fft Spectrum in the Positive and Negative Frequencies")
% xlabel("f (Hz)")
% ylabel("|fft(X)|")
% plot(time,Ft1)
%% peak
% 示例载荷时间序列 (load_time_series 为原始载荷数据)
load_time_series = Ft1;  % 原始载荷数据
% load_time_series = Ft1;  % 原始载荷数据
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

% 可视化波峰和波谷
figure;
plot(time, load_time_series);
hold on;
plot(extrema_locs_sorted, extrema_values_sorted, 'ro');
legend('原始数据', '波峰波谷');
xlabel('时间');
ylabel('载荷');
title('波峰和波谷提取');

% % 找到最高波峰和最低波谷
% [max_peak_value, max_peak_idx] = max(peaks);
% [min_valley_value, min_valley_idx] = min(valleys);
% 
% % 确定从最高波峰或最低波谷处开始
% if abs(max_peak_value) >= abs(min_valley_value)
%     start_idx = peak_locs(max_peak_idx);  % 从最高波峰开始
% else
%     start_idx = valley_locs(min_valley_idx);  % 从最低波谷开始
% end
% 
% % 找到从起始点开始的截取序列
% start_loc = find(extrema_locs_sorted >= start_idx, 1);
% extrema_locs_trimmed = extrema_locs_sorted(start_loc:end);
% extrema_values_trimmed = extrema_values_sorted(start_loc:end);
% 
% % 可视化从最高波峰/最低波谷处开始的波峰波谷序列
% figure;
% plot(time, load_time_series);
% hold on;
% plot(extrema_locs_trimmed, extrema_values_trimmed, 'ro');
% legend('原始数据', '从起始点的波峰波谷');
% xlabel('时间');
% ylabel('载荷');
% title('从起始点开始的波峰波谷');


%%  help rainflow
fs = 1;
% rainflow(Ft1,fs);
c = rainflow(extrema_values_sorted,extrema_locs_sorted);
% c2 = three_point_rainflow(load_time_series);

%goodman
theta_b = 5e7;
S = c(:,2)./(1-c(:,3)./theta_b);
%
m = 10;
N = 42565440.4361;
L
L_p = (sum(S.^m.*c(:,1)./N))^(1/m)
%
C = 9.77e70;
N = (C./S.^m);
Df
Df_p = sum(c(:,1)./N)
%%
timeL = zeros(100,1);
timeDf = zeros(100,1);
for i=2:100
    index = find(c(:,5)==i);
    if isempty(index)
       timeL(i) = timeL(i-1);
       timeDf(i) = timeDf(i-1);
    else
    %goodman
        theta_b = 5e7;
        S = c(index,2)./(1-c(index,3)./theta_b);
        %
        m = 10;
        N = 42565440.4361;
        L_p = (sum(S.^m.*c(index,1)./N))^(1/m);
        %
        C = 9.77e70;
        N = (C./S.^m);
        Df_p = sum(c(index,1)./N);
        timeL(i) = timeL(i-1)+L_p;
        timeDf(i) = timeDf(i-1)+Df_p;
    end
end

figure()
plot(1:100,timeDf)
title('Df')
figure()
plot(1:100,timeL)
title('L')