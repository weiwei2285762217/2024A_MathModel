
%% loaddata
clear;clc;
loaddata4

%%
load fitTshaft
load fitFt
load fitPitch
load fitWr
global fitTshaft fitFt fitPitch fitWr
global Pref Wind Tshaft Ft Pout Pitch Wr Wg%noise
global Pref2 Tshaft2 Ft2 Pout2 Pitch2 Wr2 Wg2%diaodu
% global Pref3 Wind3 Tshaft3 Ft3 Pout3 Pitch3 Wr3 Wg3%real
Tshaft = zeros(300,10);
Ft = zeros(300,10);
Pref2 = zeros(300,10);
Pref = Pref2;
Tshaft2 = zeros(300,10);
Ft2 = Ft;
Fs = 1; 
Fc = 0.4;
% 设计低通滤波器
d = designfilt('lowpassfir', 'FilterOrder', 1, ...
                'CutoffFrequency', Fc, ...
                'SampleRate', Fs);


for wt_num = 1:10
    Pref(:,wt_num) = Pref4(:,wt_num);
    index = find(Pref(1:end-1,wt_num)-Pref(2:end,wt_num)==0);
    for i=1:length(index)
        flag = index(i);
        if flag>10
            Pref(flag,wt_num) = mean(Pref(flag-10:flag-1,wt_num));
        else
            Pref(flag,wt_num) = mean(Pref(1:flag-1,wt_num));
        end
    end
    signal = Pref(:,wt_num);
    % 使用滤波器
    a=0;b=1-a;
    Pref(:,wt_num) = b*filter(d, signal)+a*Pref(:,wt_num);

    
    % 信号
    Wind(:,wt_num) = Wind4(:,wt_num);
    index = find(Wind(1:end-1,wt_num)-Wind(2:end,wt_num)==0);
    for i=1:length(index)
        flag = index(i);
        if flag>10
            Wind(flag,wt_num) = mean(Wind(flag-10:flag-1,wt_num));
        else
            Wind(flag,wt_num) = mean(Wind(1:flag-1,wt_num));
        end
    end
    signal = Wind(:,wt_num);
    % 使用滤波器
    a=0;b=1-a;
    Wind(:,wt_num) = b*(filter(d, signal)+a*Wind(:,wt_num));
    
    Tshaft(:,wt_num) = fitTshaft(2)+fitTshaft(3).*Pref(:,wt_num);
    Ft(:,wt_num)  = Ft4(:,wt_num);
   
    Pref2(:,wt_num) = Pref(:,wt_num);
    Tshaft2(:,wt_num)  = Tshaft(:,wt_num);
    Ft2(:,wt_num)  = Ft(:,wt_num);
end



%%

ddPref = zeros(300,100);
t=6;

tlength = 300;
for i=t:tlength

    %
    % f = coeffs(target,x(1));
    % f2 = zeros(1,100);
    % for i = 1:100
    %     f2(i) = double(f(i));
    % end

%     fun = matlabFunction(target);
%     fun_char = char(fun);
%     fun_char(4:393) = [];
%     for ichar=100:-1:1
%         charx  = 'x';
%         char1 = [charx num2str(ichar)];
%         char2 = [charx '(' num2str(ichar) ')'];
%         fun_char = strrep(fun_char, char1, char2); 
%     end
%     fun = str2func(fun_char);
    
    fun = @(x) objectfun_v4(x,i);

    
    x0 = ones(1,10)*(mean(Pref(i,:)));
%     x0 = Pref(i,:);
    % ff = feval(fun,[x0]);

    A=[];b=[];
    Aeq = ones(1,10);
    beq = sum(Pref(i,:));
    lb = ones(1,10)*(mean(Pref(i,:))-1e6);
    ub = ones(1,10)*min((mean(Pref(i,:))+1e6),5e6);
    % options = optimoptions('linprog','Algorithm','interior-point');%interior-point-legacy
    % newPref = linprog(f2,A,b,Aeq,beq,lb,ub,options);
    opts = optimoptions('fmincon','MaxFunctionEvaluations', 10000);
%     opts = optimoptions('fmincon','Algorithm','active-set'); %interior-point
    opts = optimoptions('fmincon','Algorithm','sqp');
    tic
    [newPref,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],opts);
    i
    toc
%     ddPref(i,:) = newPref;
    
    Pref2(i,:) = newPref;
    %     X1last = [Wind(i-3,:);Pref(i-3,:)];
    X1 = [Wind(i-2,:);Pref2(i-2,:)];  % 特征变量矩阵
%     X2last = [Wind(i-2,:);Pref(i-2,:)];
    X2 = [Wind(i-1,:);Pref2(i-1,:)];

%     X3last = Pref(i-1,:);
    X3 = newPref;
%     X4last = [Wind(i-1,:) ;Pitch(i-1,:) ;Wr(i-2,:)];
    X4 = [Wind(i,:)];

    Tshaftnow = fitTshaft(2)+fitTshaft(3).*X3;
    Tshaft2(i,:) = Tshaftnow;
%     Tshaftlast = fitTshaft(2)+fitTshaft(3).*X3last;

    Pitchnow = fitPitch(2)+fitPitch(3).*X1(1,:)+fitPitch(4).*X1(2,:);
%     Pitchlast = fitPitch(2)+fitPitch(3).*X1last(1,:)+fitPitch(4).*X1last(2,:);
    Wrnow = fitWr(2)+fitWr(3).*X2(1,:)+fitWr(4).*X2(2,:);
%     Wrlast = fitWr(2)+fitWr(3).*X2last(1,:)+fitWr(4).*X2last(2,:);

    Ftnow = fitFt(2)+fitFt(3).*X4(1,:)+fitFt(4).*Pitchnow+fitFt(5).*Wrnow;
    Ft2(i,:) = Ftnow;
    % [x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
%     toc
    
end
save('Pref2_v4.mat','Pref2')
save('Tshaft2_v4.mat','Tshaft2')
save('Ft2_v4.mat','Ft2')
%'interior-point' (default)'trust-region-reflective''sqp''sqp-legacy' (optimoptions only)'active-set'
%%
figure()
load('Pref2_v4.mat','Pref2')
load('Tshaft2_v4.mat','Tshaft2')
load('Ft2_v4.mat','Ft2')
Time = 1:300;
n=1;tlength = 2000;
% chaPref = (Pref(:,n)-Pref2(:,n));
plot(Time,Pref(:,n))
hold on
plot(Time,Pref2(:,n))
plot(Time,min((mean(Pref,2)+1e6),5e6))
plot(Time,(mean(Pref,2)-1e6))
plot(Time,(mean(Pref,2)))
xlim([0 100])
% ylim([-1 1])

table = zeros(2,6);
table(1,1) = std(Pref(:,n));
table(2,1) = mean(Pref(:,n));

table(1,2) = std(Pref2(:,n));
table(2,2) = mean(Pref2(:,n));
legend('干扰下原始','调度','上','下','中')
% mean(Wind(:,n))
% mean(Pref2(:,n))
%%
figure()
Time = 1:300;
plot(Time,Tshaft(:,n))
hold on
plot(Time,Tshaft2(:,n))
title('Tshaft')
legend('原始','新')

table(1,3) = std(Tshaft(:,n));
table(2,3) = mean(Tshaft(:,n));
table(1,4) = std(Tshaft2(:,n));
table(2,4) = mean(Tshaft2(:,n));

figure()
plot(Time,Ft(:,n))
hold on
plot(Time,Ft2(:,n))
title('Ft')
legend('原始','新')

table(1,5) = std(Ft(:,n));
table(2,5) = mean(Ft(:,n));
table(1,6) = std(Ft2(:,n));
table(2,6) = mean(Ft2(:,n));
table


%%
t=10;
% chaPref = Pref(:,n)-ddPref(:,n);
plot(1:10,Pref(t,:))
hold on
plot(1:10,Pref2(t,:))



% xlim([0 tlength])