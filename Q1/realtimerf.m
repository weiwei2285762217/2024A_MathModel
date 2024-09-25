%%
allL = zeros(1,100);
allDf = zeros(1,100);
TimeL = zeros(101,100);
TimeDf = zeros(101,100);
prev_stress = [0];
tic
for i=1:100
    Ft1 = Ft(1:100,i);
    load_time_series = Ft1;  % 原始载荷数据
    j=3;
    prev_stress = load_time_series(j-2:j-1);
    for j = 3:1:100
        
        [c,prev_stress] = realtime_three_point_rainflow(load_time_series(j),prev_stress);
        if isempty(c)
           c=[0 1 1]; 
        end
        %goodman
        theta_b = 5e7;
        S = c(:,2)./(1-c(:,3)./theta_b);
        %
        m = 10;
        N = 42565440.4361;
    %     L
        TimeL(j+1,i) = (sum(S.^m.*c(:,1)./N))^(1/m) + TimeL(j,i);
        %
        C = 9.77e70;
        N = (C./S.^m);
    %     Df
        TimeDf(j+1,i) = sum(c(:,1)./N) + TimeL(j,i);
    end
    
end
L = Ft(101,:);
Df = Ft(102,:);
toc
figure()

plot(1:100,TimeL(101,:))
hold on
plot(1:100,L)
title('L')
legend('pre','real')
figure()

plot(1:100,TimeDf(101,:))
hold on
plot(1:100,Df)
legend('pre','real')
title('Df')

figure()
plot(1:100,TimeL(1:100,1))
title('TimeL1')
figure()
plot(1:100,TimeDf(1:100,1))
title('TimeDf')