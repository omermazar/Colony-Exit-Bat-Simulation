function [] =Test()

MinDist = 0.05;
t= 0:0.01:5;
x0Prey = 5;
y0Prey = 5;
V0Prey= 0.1;
xPrey= [x0Prey, zeros(1,length(t)-1)];
yPrey= [y0Prey, zeros(1,length(t)-1)];
Teta0Prey= -0;
DTetaPrey= zeros(1,length(t)); % randn(1,length(t))*pi/36;
TetaPrey = [Teta0Prey,  zeros(1,length(t)-1)];

V0Hunter= 0.4;
x0Hunter = 0;
y0Hunter = 0;
xHunterPT = [x0Hunter, zeros(1,length(t)-1)];
yHunterPT = [y0Hunter, zeros(1,length(t)-1)];
xHunterPN = xHunterPT;
yHunterPN = yHunterPT;

Lamda0 =   atan((y0Prey-y0Hunter)/(x0Prey-x0Hunter));
Distance0 = sqrt((y0Prey-y0Hunter)^2+(x0Prey-x0Hunter)^2);
Delta0 = pi/36; % asin(V0Prey/V0Hunter*sin(TetaPrey));
Lamda = [Lamda0, zeros(1,length(t)-1)];
Delta = [Delta0, zeros(1,length(t)-1)];
DistancePT = [Distance0, -1*ones(1,length(t)-1)];
DistancePN = [Distance0, -1*ones(1,length(t)-1)];
MinDistanceAchievedPT =   Distance0;
MinDistanceAchievedPN = Distance0;

TetaHunterPN = [Lamda0+ Delta0, zeros(1,length(t)-1)]; % Parralel Navigation
TetaHunterPT = [Lamda0, zeros(1,length(t)-1)]; % Pure TRcking

CatchFlagPT = 0; MissFlagPT=0;
xCatchPT = [];
yCatchPT = [];
TimeOfCatchPT=[];
TimeOfMissPT =[];xMissPT=[]; yMissPT=[];

CatchFlagPN = 0; MissFlagPN=0;
xCatchPN = [];
yCatchPN = [];
TimeOfCatchPN=[];
TimeOfMissPN =[];xMissPN=[]; yMissPN=[];


VHunter = V0Hunter*ones(1,length(t));

DeltaPT=0;
RPN(1) = Distance0;
dRPN(1) =0;
%%%%% Pure Tracking
for n= 1:length(t)-1
    DistancePT(n) = sqrt( (yPrey(n)-yHunterPT(n))^2 + (xPrey(n)-xHunterPT(n))^2 );
%     DistancePN(n) = sqrt( (yPrey(n)-yHunterPN(n))^2 + (xPrey(n)-xHunterPN(n))^2 );
    if DistancePT(n)<= MinDist % Catch
        CatchFlagPT =1;
        TimeOfCatchPT = [TimeOfCatchPT, n];
        xCatchPT = [xCatchPT, xHunterPT(n)];
        yCatchPT = [yCatchPT, yHunterPT(n)];
%         break
    end % Distace
    % MISS
    MinDistanceAchievedPT = min(MinDistanceAchievedPT , DistancePT(n));
    DDistPT =DistancePT(n)- DistancePT(max(1,n-1));
    
    xPrey(n+1) = xPrey(n)+ V0Prey*cos(TetaPrey(n));
    yPrey(n+1) = yPrey(n)+ V0Prey*sin(TetaPrey(n));
    TetaPrey(n+1) = TetaPrey(n)+ DTetaPrey(n);
    
    if DDistPT >0 && DistancePT(n) >= 2*MinDistanceAchievedPT % The distnce encreeasing
        MissFlagPT = 1;
        TimeOfMissPT = [TimeOfMissPT, n];
        xMissPT = [xMissPT, xHunterPT(n)];
        yMissPT = [yMissPT, yHunterPT(n)];
%         break
    end % DDIST
    
    
    LamdaPT(n+1) = atan((yPrey(n)-yHunterPT(n))/(xPrey(n)-xHunterPT(n))); % the Spatial angle betwwen prey and Hunter
    TetaHunterPT(n+1) = LamdaPT(n+1);
    xHunterPT(n+1) = xHunterPT(n)+ VHunter(n)*cos(TetaHunterPT(n)-DeltaPT);
    yHunterPT(n+1) = yHunterPT(n)+ VHunter(n)*sin(TetaHunterPT(n)-DeltaPT);
    
    DistancePN(n) = sqrt( (yPrey(n)-yHunterPN(n))^2 + (xPrey(n)-xHunterPN(n))^2 );
    if DistancePN(n)<= MinDist % Catch
        CatchFlagPN =1;
        TimeOfCatchPN = [TimeOfCatchPN, n];
        xCatchPN = [xCatchPN, xHunterPN(n)];
        yCatchPN = [yCatchPN, yHunterPN(n)];
        break
    end % Distace
      % MISS
    
    MinDistanceAchievedPN = min(MinDistanceAchievedPN , DistancePN(n));
    DDistPN =DistancePN(n)- DistancePN(max(1,n-1));
    if DDistPN >0 && DistancePN(n) >= 2*MinDistanceAchievedPN % The distnce encreeasing
        MissFlagPN = 1;
        TimeOfMissPN = [TimeOfMissPN, n];
        xMissPN = [xMissPN, xHunterPN(n)];
        yMissPN = [yMissPN, yHunterPN(n)];
        break
    end % DDIST
% % % % % % % % % % %     Teta(n) = TetaPrey(n) - Lamda0;
% % % % % % % % % % %     drPN(n) = V0Prey*cos(Teta(n))-V0Hunter*cos(Delta(n));
% % % % % % % % % % %     rPN(n+1) = rPN(n)+drPN(n);
% % % % % % % % % % %     dLamda(n) = 1./rPN(n)*(V0Prey*sin(Teta(n))-V0Hunter*sin(Delta(n)));
% % % % % % % % % % %     Delta(n+1) = Delta(n)+ dLamda(n);
% % % % % % % % % % % %     dDelta(n) =  -1*dLamda(n);    
% % % % % % % % % % %     TetaHunterPN(n+1) = Lamda0+ Delta(n+1); %N=1
% % % % % % % % % % % %     Delta(n+1) = Delta(n)+ dDelta(n);
% % % % % % % % % % % %     dTetaHunter(n) = -2.25*V0Hunter / Distance0 * sin(Lamda0) * (DistancePN(n)/Distance0)^1;
% % % % % % % % % % % %     TetaHunterPN(n+1) = TetaHunterPN(n)+dTetaHunterPN(n);
% % % % % % % % % % % %     LamdaPN(n+1) = atan((yPrey(n)-yHunterPN(n))/(xPrey(n)-xHunterPN(n)));
% % % % % % % % % % % %     dLamda(n) = LamdaPN(n+1)- LamdaPN(n);
% % % %     Teta(n) = TetaPrey(n) - Lamda0;
% % % %     dRPN(n) = V0Prey*cos(Teta(n)) - V0Hunter*cos(Delta(n));
% % % %     RPN(n+1) = RPN(n)+ dRPN(n);
% % % %     dLamda(n) = 1/DistancePN(n) * (V0Prey*sin(Teta(n)) - V0Hunter*sin(Delta(n)));
% % % %     Delta(n+1) = Delta(n)+dLamda(n); 
% % % %     TetaHunterPN(n+1) = TetaHunterPN(n)+ dLamda(n);
    
    Lamda(n) =   atan((yPrey(n)-yHunterPN(n))/(xPrey(n)-xHunterPN(n)));
    dDelta(n) =  Lamda0- Lamda(n);
    TetaHunterPN(n+1) = TetaHunterPN(n) - dDelta(n);
    xHunterPN(n+1) = xHunterPN(n)+ VHunter(n)*cos(TetaHunterPN(n+1));
    yHunterPN(n+1) = yHunterPN(n)+ VHunter(n)*sin(TetaHunterPN(n+1));
end % nn

% xPrey = xPrey(IndPT);
% yPrey = yPrey(IndPT);


IndPT = find(DistancePT >= 0);
DistancePT = DistancePT(IndPT);
xHunterPT = xHunterPT(IndPT);
yHunterPT = yHunterPT(IndPT);
TetaHunterPT = TetaHunterPT(IndPT);

IndPN = find(DistancePN >= 0);
DistancePN = DistancePN(IndPN);
xHunterPN = xHunterPN(IndPN);
yHunterPN = yHunterPN(IndPN);
TetaHunterPN = TetaHunterPN(IndPT);



figure
hold on
plot(xPrey(1:n),yPrey(1:n))
plot(xHunterPT,yHunterPT,'r')
if CatchFlagPT
    plot(xCatchPT,yCatchPT,'r*')
end  
if MissFlagPT
    plot(xMissPT,yMissPT,'ro')
end

plot(xHunterPN,yHunterPN,'k')
if CatchFlagPN
    plot(xCatchPN,yCatchPN,'k*')
end  
if MissFlagPN
    plot(xMissPN,yMissPN,'ko')
end


end
