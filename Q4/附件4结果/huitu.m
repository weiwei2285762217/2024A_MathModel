
%% Tshaft - -realtime
load L_NT1.mat 
load D_FT1.mat 
load L_NT2.mat 
load D_FT2.mat
%%
figure()
plot(1:10,L_NT1(300,:))
hold on
plot(1:10,L_NT2(300,:))
legend('优化前','优化后')

sum(L_NT1(300,:))
sum(L_NT2(300,:))

figure()
plot(1:10,log10(D_FT1(300,:)))
hold on
plot(1:10,log10(D_FT2(300,:)))
legend('优化前','优化后')
% ylim([0 2e-17])
sum(D_FT1(300,:))
sum(D_FT2(300,:))
%%
figure()
plot(1:10,D_FT2(:,50))

%% Ft - -realtime
load L_NF1.mat 
load D_FF1.mat 
load L_NF2.mat 
load D_FF2.mat
%%
figure()
plot(1:10,L_NF1(300,:))
hold on
plot(1:10,L_NF2(300,:))
legend('优化前','优化后')
sum(L_NF1(300,:))
sum(L_NF2(300,:))
figure()
plot(1:10,log10(D_FF1(300,:)))
hold on
plot(1:10,log10(D_FF2(300,:)))
legend('优化前','优化后')
% ylim([0 2e-17])
sum(D_FF1(300,:))
sum(D_FF2(300,:))
%%
figure()
plot(1:10,D_FF2(:,50))
%% 优化前后风机间参考功率方差值
load('Pref2_v2.mat')
Pref3_std = std(Pref3');
Pref2_std = std(Pref2');
Pref_std = std(Pref');
figure()

plot(1:2000,Pref3_std.^2)
hold on
plot(1:2000,Pref_std.^2)
plot(1:2000,Pref2_std.^2)
legend('原始','干扰','优化')

