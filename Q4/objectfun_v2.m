function target = objectfun(x,i)
    %     i=t+1;%i-1刚好当前
    global Pref Wind Tshaft Ft Pout Pitch Wr Wg
    global fitTshaft fitFt fitPitch fitWr
    global Pref2 Tshaft2 Ft2 Pout2 Pitch2 Wr2 Wg2
    Tshaftlast = Tshaft2(i-1:i-1,:);
    Ftlast = Ft2(i-1:i-1,:);
    %now
%     X1last = [Wind(i-3,:);Pref(i-3,:)];
    X1 = [Wind(i-2,:);Pref2(i-2,:)];  % 特征变量矩阵
%     X2last = [Wind(i-2,:);Pref(i-2,:)];
    X2 = [Wind(i-1,:);Pref2(i-1,:)];
%     X3last = Pref(i-1,:);
    X3 = x;
%     X4last = [Wind(i-1,:) ;Pitch(i-1,:) ;Wr(i-2,:)];
    X4 = [Wind(i,:) ;Pitch(i,:) ;Wr(i-1,:)];
    Tshaftnow = fitTshaft(2)+fitTshaft(3).*X3;
%     Tshaftlast = fitTshaft(2)+fitTshaft(3).*X3last;
    Pitchnow = fitPitch(2)+fitPitch(3).*X1(1,:)+fitPitch(4).*X1(2,:);
%     Pitchlast = fitPitch(2)+fitPitch(3).*X1last(1,:)+fitPitch(4).*X1last(2,:);
    Wrnow = fitWr(2)+fitWr(3).*X2(1,:)+fitWr(4).*X2(2,:);
%     Wrlast = fitWr(2)+fitWr(3).*X2last(1,:)+fitWr(4).*X2last(2,:);
    Ftnow = fitFt(2)+fitFt(3).*X4(1,:)+fitFt(4).*Pitchnow+fitFt(5).*Wrnow;
%     Ftlast = fitFt(2)+fitFt(3).*X4last(1,:)+fitFt(4).*Pitchlast+fitFt(5).*Wrlast;
    % Ftpre = fitFt(2)+fitFt(3).*X4(1,:)+fitFt(4).*X4(2,:)+fitFt(5).*X4(3,:);
    % now+1
    %     X1last = [Wind(i-3,:);Pref(i-3,:)];
    X1 = [Wind(i-1,:);Pref2(i-1,:)];  % 特征变量矩阵
%     X2last = [Wind(i-2,:);Pref(i-2,:)];
    X2 = [Wind(i,:);x];
    X3 = x;
    X4 = [Wind(i,:) ;Pitch(i,:) ;Wr(i,:)];
    Tshaftnow_2 = fitTshaft(2)+fitTshaft(3).*X3;
    Pitchnow_2 = fitPitch(2)+fitPitch(3).*X1(1,:)+fitPitch(4).*X1(2,:);
    Wrnow_2 = fitWr(2)+fitWr(3).*X2(1,:)+fitWr(4).*X2(2,:);
    Ftnow_2 = fitFt(2)+fitFt(3).*X4(1,:)+fitFt(4).*Pitchnow_2+fitFt(5).*Wrnow_2;
    % now+2
    %     X1last = [Wind(i-3,:);Pref(i-3,:)];
    X1 = [Wind(i,:);x];  % 特征变量矩阵
%     X2last = [Wind(i-2,:);Pref(i-2,:)];
    X2 = [Wind(i,:);x];
%     X3last = Pref(i-1,:);
    X3 = x;
%     X4last = [Wind(i-1,:) ;Pitch(i-1,:) ;Wr(i-2,:)];
    X4 = [Wind(i,:) ;Pitch(i,:) ;Wr(i,:)];
    Tshaftnow_3 = fitTshaft(2)+fitTshaft(3).*X3;
%     Tshaftlast = fitTshaft(2)+fitTshaft(3).*X3last;
    Pitchnow_3 = fitPitch(2)+fitPitch(3).*X1(1,:)+fitPitch(4).*X1(2,:);
%     Pitchlast = fitPitch(2)+fitPitch(3).*X1last(1,:)+fitPitch(4).*X1last(2,:);
    Wrnow_3 = fitWr(2)+fitWr(3).*X2(1,:)+fitWr(4).*X2(2,:);
%     Wrlast = fitWr(2)+fitWr(3).*X2last(1,:)+fitWr(4).*X2last(2,:);
    Ftnow_3 = fitFt(2)+fitFt(3).*X4(1,:)+fitFt(4).*Pitchnow_3+fitFt(5).*Wrnow_3;
%     Ftlast = fitFt(2)+fitFt(3).*X4last(1,:)+fitFt(4).*Pitchlast+fitFt(5).*Wrlast;
    % Ftpre = fitFt(2)+fitFt(3).*X4(1,:)+fitFt(4).*X4(2,:)+fitFt(5).*X4(3,:);
    
    
    a=1;b=1e2;c=1e2;
%     Tshaft_target = [Tshaftlast-Tshaftnow]./mean(Tshaftlast);
    Tshaft_target = [Tshaftnow_3-Tshaftlast ; Tshaftnow_2-Tshaftlast ; Tshaftnow-Tshaftlast ];%./mean(Tshaftlast);
    Tshaft_target = mean(Tshaft_target);
    Ft_target = [Ftnow_3-Ftlast ; Ftnow_2-Ftlast ; Ftnow-Ftlast];
    Ft_target = mean(Ft_target);
    Prefcha = abs(x-Pref(i-1,:));
    % target = sum(a*Tshaftpre+b*Ftpre);
    target = sum(a*(Tshaft_target).^2+b*(Ft_target).^2+c*Prefcha,'all');
end
