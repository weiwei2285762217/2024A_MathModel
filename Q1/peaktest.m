
%% Tq
opts = spreadsheetImportOptions("NumVariables", 100);

% 指定工作表和范围
opts.Sheet = "Sheet1";
opts.DataRange = "B2:CW103";

% 指定列名称和类型
opts.VariableNames = ["WT1", "WT2", "WT3", "WT4", "WT5", "WT6", "WT7", "WT8", "WT9", "WT10", "WT11", "WT12", "WT13", "WT14", "WT15", "WT16", "WT17", "WT18", "WT19", "WT20", "WT21", "WT22", "WT23", "WT24", "WT25", "WT26", "WT27", "WT28", "WT29", "WT30", "WT31", "WT32", "WT33", "WT34", "WT35", "WT36", "WT37", "WT38", "WT39", "WT40", "WT41", "WT42", "WT43", "WT44", "WT45", "WT46", "WT47", "WT48", "WT49", "WT50", "WT51", "WT52", "WT53", "WT54", "WT55", "WT56", "WT57", "WT58", "WT59", "WT60", "WT61", "WT62", "WT63", "WT64", "WT65", "WT66", "WT67", "WT68", "WT69", "WT70", "WT71", "WT72", "WT73", "WT74", "WT75", "WT76", "WT77", "WT78", "WT79", "WT80", "WT81", "WT82", "WT83", "WT84", "WT85", "WT86", "WT87", "WT88", "WT89", "WT90", "WT91", "WT92", "WT93", "WT94", "WT95", "WT96", "WT97", "WT98", "WT99", "WT100"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% 导入数据
Tq = readtable("H:\DATA\A\Q1\Tq.xlsx", opts, "UseExcel", false);

Tq = table2array(Tq);

clear opts
%% Ft
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
n=80;
time = 1:1:100;
Ft1 = Ft(1:100,n);
L = Ft(101,n);
Df = Ft(102,n);

%% 频域
mainFs = zeros(1,100);
for i=1:100
    Fs = 1;
    n=i;
    X = Ft(1:100,n)-mean(Ft(1:100,n));%0点偏移
    % X = highpass(Ft1, 0.1, Fs);  % fs 是采样频率
    len = length(X);
    Y = fft(X);
    P2 = abs(Y/len);
    P1 = P2(1:len/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    [~, maxindex] = max(P1);
    f = Fs/len*(0:len/2);
    mainFs(1,i) = f(maxindex);
%     , Fs/len*(-len/2:len/2-1)
end
plot(1:100,mainFs)
%% 峰谷均值分析
meanpeak = zeros(1,100);
meanvalley = meanpeak;
meanload = meanpeak;
for i=1:100
    n=i;
    load_time_series = Ft(1:100,n);  % 原始载荷数据
    time = 1:length(load_time_series);  % 时间序列

    % 提取波峰
    [peaks, peak_locs] = findpeaks(load_time_series, time);

    % 提取波谷 (通过对数据取反找波峰)
    [valleys, valley_locs] = findpeaks(-load_time_series, time);
    valleys = -valleys;  % 将波谷恢复为原始值
    meanpeak(1,i) = mean(peaks);
    meanvalley(1,i) = mean(valleys);
    meanload(1,i) = mean(load_time_series); 
end
plot(1:100,meanpeak,1:100,meanvalley,1:100,meanload)
legend('peak','valley','load')
%% 峰谷分析
allpeaks = [];
allvalleys = [];
for i=1:100
    n=i;
    load_time_series = Ft(1:100,n);  % 原始载荷数据
    time = 1:length(load_time_series);  % 时间序列

    % 提取波峰
    [peaks, peak_locs] = findpeaks(load_time_series, time);

    % 提取波谷 (通过对数据取反找波峰)
    [valleys, valley_locs] = findpeaks(-load_time_series, time);
    valleys = -valleys;  % 将波谷恢复为原始值
    allpeaks = [allpeaks;peaks];
    allvalleys = [allvalleys;valleys];
end
figure(1)
plot(1:length(allpeaks),allpeaks)
hold on
plot(1:length(allvalleys),allvalleys)
legend('peak','valley')
