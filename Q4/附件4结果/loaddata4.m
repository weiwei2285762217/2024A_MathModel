%% 设置导入选项并导入数据
clc;clear
opts = spreadsheetImportOptions("NumVariables", 10);
opts.Sheet = "Ft";
opts.DataRange = "A1:J300";
opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
Ft4 = readtable("H:\DATA\A\Q4\附件4-噪声和延迟作用下的采集数据.xlsx", opts, "UseExcel", false);
Ft4 = table2array(Ft4);
clear opts
%% 设置导入选项并导入数据
opts = spreadsheetImportOptions("NumVariables", 10);
opts.Sheet = "Pref";
opts.DataRange = "A1:J300";
opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
Pref4 = readtable("H:\DATA\A\Q4\附件4-噪声和延迟作用下的采集数据.xlsx", opts, "UseExcel", false);
Pref4 = table2array(Pref4);
clear opts
%% 设置导入选项并导入数据
opts = spreadsheetImportOptions("NumVariables", 10);
opts.Sheet = "Vwin";
opts.DataRange = "A1:J300";
opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
Wind4 = readtable("H:\DATA\A\Q4\附件4-噪声和延迟作用下的采集数据.xlsx", opts, "UseExcel", false);
Wind4 = table2array(Wind4);
clear opts
