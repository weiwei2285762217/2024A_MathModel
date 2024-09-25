%% loaddata
clear;clc;
load('附件2-风电机组采集数据.mat','data_TS_WF');
WT1 = data_TS_WF.WF_1.WT; 
WT2 = data_TS_WF.WF_2.WT;
%%
WT=WT1;
for wt_num = 1:100
    tempwt = WT{1,wt_num};
    Time = tempwt.time;
    Dt = Time(2)-Time(1);
    Pref(:,wt_num) = tempwt.inputs(:,1);
    Wind(:,wt_num)  = tempwt.inputs(:,2);
    Tshaft(:,wt_num)  = tempwt.outputs(:,1);
    Ft(:,wt_num)  = tempwt.outputs(:,2);
    Pout(:,wt_num)  = tempwt.outputs(:,3);
    Pitch(:,wt_num)  = tempwt.states(:,1);
    Wr(:,wt_num)  = tempwt.states(:,2);
    Wg(:,wt_num)  = tempwt.states(:,3);
%     WTdata = [Pref Wind Tshaft Ft Pout Pitch Wr Wg];
end
%% test
allPref = zeros(2000,100);
for i=1:100
    tempwt = WT{1,i};
    Pref = tempwt.inputs(:,1);
    allPref(:,i) = Pref;
end
allPref_mean = mean(allPref,2);
figure()
plot(1:2000,allPref_mean)
title('Perfmean')

allPref_sum = sum(allPref,2);
figure()
plot(1:2000,allPref_sum)
title('Perfsum')
%%
tic
ddPref = zeros(2000,100);
syms x [1 100]
t=5;
load fitTshaft
load fitFt
load fitPitch
load fitWr
tlength = 50;
for i=t:tlength
%     i=t+1;%i-1刚好当前
    X1last = [Wind(i-3,:);Pref(i-3,:)];
    X1 = [Wind(i-2,:);Pref(i-2,:)];  % 特征变量矩阵
    X2last = [Wind(i-2,:);Pref(i-2,:)];
    X2 = [Wind(i-1,:);Pref(i-1,:)];

    X3last = Pref(i-1,:);
    X3 = x;
    X4last = [Wind(i-1,:) ;Pitch(i-1,:) ;Wr(i-2,:)];
    X4 = [Wind(i,:) ;Pitch(i,:) ;Wr(i-1,:)];

    Tshaftnow = fitTshaft(2)+fitTshaft(3).*X3;
    Tshaftlast = fitTshaft(2)+fitTshaft(3).*X3last;

    Pitchnow = fitPitch(2)+fitPitch(3).*X1(1,:)+fitPitch(4).*X1(2,:);
    Pitchlast = fitPitch(2)+fitPitch(3).*X1last(1,:)+fitPitch(4).*X1last(2,:);
    Wrnow = fitWr(2)+fitWr(3).*X2(1,:)+fitWr(4).*X2(2,:);
    Wrlast = fitWr(2)+fitWr(3).*X2last(1,:)+fitWr(4).*X2last(2,:);

    Ftnow = fitFt(2)+fitFt(3).*X4(1,:)+fitFt(4).*Pitchnow+fitFt(5).*Wrnow;
    Ftlast = fitFt(2)+fitFt(3).*X4last(1,:)+fitFt(4).*Pitchlast+fitFt(5).*Wrlast;
    % Ftpre = fitFt(2)+fitFt(3).*X4(1,:)+fitFt(4).*X4(2,:)+fitFt(5).*X4(3,:);

    a=0.01;b=1-a;

    % target = sum(a*Tshaftpre+b*Ftpre);
    target = sum(a*abs(Tshaftnow-Tshaftlast)+b*abs(Ftnow-Ftlast));
    %
    % f = coeffs(target,x(1));
    % f2 = zeros(1,100);
    % for i = 1:100
    %     f2(i) = double(f(i));
    % end

    fun = matlabFunction(target);
    fun_char = char(fun);
    fun_char(4:393) = [];
    for ichar=100:-1:1
        charx  = 'x';
        char1 = [charx num2str(ichar)];
        char2 = [charx '(' num2str(ichar) ')'];
        fun_char = strrep(fun_char, char1, char2); 
    end
    fun = str2func(fun_char);
    x0 = Pref(i,:);
    % ff = feval(fun,[x0]);

    A=[];b=[];
    Aeq = ones(1,100);
    beq = sum(Pref(i,:));
    lb = ones(1,100)*(mean(Pref(i,:))-1e6);
    ub = ones(1,100)*min((mean(Pref(i,:))+1e6),5e6);
    % options = optimoptions('linprog','Algorithm','interior-point');%interior-point-legacy
    % newPref = linprog(f2,A,b,Aeq,beq,lb,ub,options);
    opts = optimoptions('fmincon','MaxFunctionEvaluations', 10000);
%     opts = optimoptions('fmincon','Algorithm','active-set'); %interior-point
    opts = optimoptions('fmincon','Algorithm','interior-point');
%     tic
    [newPref,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],opts);
    i
    ddPref(i,:) = newPref;
    % [x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
%     toc
    
end
toc
%'interior-point' (default)'trust-region-reflective''sqp''sqp-legacy' (optimoptions only)'active-set'
%%
chaPref = Pref(:,1)-ddPref(:,1);
plot(Time,Pref(:,1))
hold on
plot(Time,ddPref(:,1))
plot(Time,max((mean(Pref,2)+1e6),5e6))
plot(Time,(mean(Pref,2)-1e6))
xlim([0 tlength])
% ylim([-1 1])

legend('原始','调度','上','下')