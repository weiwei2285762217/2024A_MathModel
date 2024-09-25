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
%%

Ft_std = std(Ft(1:100,:),0,1);
Tshaft_std = std(Tshaft(1:100,:),0,1);

Ft_mean = mean(Ft(1:100,:),1);
Tshaft_mean = mean(Tshaft(1:100,:),1);

Tshaft_L = Tshaft(101,:);
Ft_L = Ft(101,:);

Tshaft_D = Tshaft(102,:);
Ft_D = Ft(102,:);
%%
corr(Ft_std',Ft_L')
corr(Ft_std',Ft_D')
corr(Tshaft_std',Tshaft_L')
corr(Tshaft_std',Tshaft_D')

%%
corr(Ft_mean',Ft_L')
corr(Ft_mean',Ft_D')
corr(Tshaft_mean',Tshaft_L')
corr(Tshaft_mean',Tshaft_D')
%用均值和L能相关

%%
corr(Ft_std',Ft_L')
corr(Ft_std',Tshaft_L')

corr(Tshaft_std',Ft_L')
corr(Tshaft_std',Tshaft_L')
%%
corr(Ft_mean',Ft_L')
corr(Ft_mean',Tshaft_L')

corr(Tshaft_mean',Ft_L')
corr(Tshaft_mean',Tshaft_L')

corr(Ft_mean',Ft_D')
corr(Ft_mean',Tshaft_D')

corr(Tshaft_mean',Ft_D')
corr(Tshaft_mean',Tshaft_D')