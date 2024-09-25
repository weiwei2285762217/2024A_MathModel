%% loaddata
clear;clc;
load('附件3-噪声和延迟作用下的采集数据.mat','data_TS_WF');
WT1 = data_TS_WF.WF_1.WT; 
WT2 = data_TS_WF.WF_2.WT;
load('附件2-风电机组采集数据.mat')
WT3 = data_TS_WF.WF_1.WT;

%%
global Pref Wind Tshaft Ft Pout Pitch Wr Wg
global Pref2 Tshaft2 Ft2 Pout2 Pitch2 Wr2 Wg2
global Pref3 Tshaft3 Ft3 Pout3 Pitch3 Wr3 Wg3
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
    Pref2(:,wt_num) = tempwt.inputs(:,1);
    Tshaft2(:,wt_num)  = tempwt.outputs(:,1);
    Ft2(:,wt_num)  = tempwt.outputs(:,2);
    Pout2(:,wt_num)  = tempwt.outputs(:,3);
    Pitch2(:,wt_num)  = tempwt.states(:,1);
    Wr2(:,wt_num)  = tempwt.states(:,2);
    Wg2(:,wt_num)  = tempwt.states(:,3);
    %
    tempwt = WT3{1,wt_num};
    Pref3(:,wt_num) = tempwt.inputs(:,1);
    Wind3(:,wt_num)  = tempwt.inputs(:,2);
    Tshaft3(:,wt_num)  = tempwt.outputs(:,1);
    Ft3(:,wt_num)  = tempwt.outputs(:,2);
    Pout3(:,wt_num)  = tempwt.outputs(:,3);
    Pitch3(:,wt_num)  = tempwt.states(:,1);
    Wr3(:,wt_num)  = tempwt.states(:,2);
    Wg3(:,wt_num)  = tempwt.states(:,3);
end

%%

ddPref = zeros(2000,100);
t=6;
load fitTshaft
load fitFt
load fitPitch
load fitWr
tlength = 2000;
global fitTshaft fitFt fitPitch fitWr
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
    
    fun = @(x) objectfun(x,i);

    
    x0 = ones(1,100)*(mean(Pref(i,:)));
%     x0 = Pref(i,:);
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
    opts = optimoptions('fmincon','Algorithm','sqp');
    tic
    [newPref,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],opts);
    i
    toc
    ddPref(i,:) = newPref;
    %加噪声与延迟
%     noise1 = normrnd(0, 0.25);
%     Pref2(i,:) = newPref.*(1+noise1);
%     if i>40 && i<50
%         newPref(1,1) = Pref2(i,1);
%     end
    %     X1last = [Wind(i-3,:);Pref(i-3,:)];
    X1 = [Wind(i-2,:);Pref2(i-2,:)];  % 特征变量矩阵
%     X2last = [Wind(i-2,:);Pref(i-2,:)];
    X2 = [Wind(i-1,:);Pref2(i-1,:)];

%     X3last = Pref(i-1,:);
    X3 = newPref;
%     X4last = [Wind(i-1,:) ;Pitch(i-1,:) ;Wr(i-2,:)];
    X4 = [Wind(i,:) ;Pitch(i,:) ;Wr(i-1,:)];

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
save('Pref2.mat','Pref2')
save('Tshaft2.mat','Tshaft2')
save('Ft2.mat','Ft2')
%'interior-point' (default)'trust-region-reflective''sqp''sqp-legacy' (optimoptions only)'active-set'
%%

load('Pref2.mat','Pref2')
load('Tshaft2.mat','Tshaft2')
load('Ft2.mat','Ft2')
n=2;tlength = 2000;
chaPref = (Pref(:,n)-Pref2(:,n));
plot(Time,Pref3(:,n))
hold on
plot(Time,Pref2(:,n))
plot(Time,min((mean(Pref3,2)+1e6),5e6))
plot(Time,(mean(Pref3,2)-1e6))
plot(Time,(mean(Pref3,2)))
xlim([0 100])
% ylim([-1 1])

table = zeros(2,6);
table(1,1) = std(Pref3(:,n));
table(2,1) = mean(Pref3(:,n));

table(1,2) = std(Pref2(:,n));
table(2,2) = mean(Pref2(:,n));
legend('原始','调度','上','下','中')
mean(Wind(:,n))
mean(Pref2(:,n))
%%
figure()
Time = 1:2000;
plot(Time,Tshaft3(:,n))
hold on
plot(Time,Tshaft2(:,n))
title('Tshaft')
legend('原始','新')

table(1,3) = std(Tshaft3(:,n));
table(2,3) = mean(Tshaft3(:,n));
table(1,4) = std(Tshaft2(:,n));
table(2,4) = mean(Tshaft2(:,n));

figure()
plot(Time,Ft3(:,n))
hold on
plot(Time,Ft2(:,n))
title('Ft')
legend('原始','新')

table(1,5) = std(Ft3(:,n));
table(2,5) = mean(Ft3(:,n));
table(1,6) = std(Ft2(:,n));
table(2,6) = mean(Ft2(:,n));
table


%%
t=10;
% chaPref = Pref(:,n)-ddPref(:,n);
plot(1:100,Pref3(t,:))
hold on
plot(1:100,Pref2(t,:))



% xlim([0 tlength])