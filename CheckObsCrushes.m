function [isCrushed, xRecover, yRecover, tetaRecover] = CheckObsCrushes(BAT, Terrain, minDistTH, currPulseNum, nStart)

% this function check whther the bat ha crushed into an obstacle during the
% last manuever and reuens the new position of the bat
%  Terrain = [946 946 956 956; 520 540 540 520  ; ...
%             948 950 955 955;  523 530 530 523 ];

%%% Changed 07Jan23

nEnd = min(nStart + BAT.TransmittedPulsesStruct(currPulseNum).IPItoNextPulse-1, numel(BAT.xBati));
% % % nStart = nTime - BAT.TransmittedPulsesStruct(currPulseNum).IPItoNextPulse +1;
sz = size(Terrain,1);

% init ouptputs
isCrushed   = false; 
xRecover    = []; % the stating pot of thr Call
yRecover    = [];
tetaRecover = [];

%%%% Current Position pf the bat
xB = BAT.xBati(nEnd);
yB = BAT.yBati(nEnd);
%%%%% the Manuever During last Call Intrval
xMan = BAT.xBati(nStart:nEnd);
yMan = BAT.yBati(nStart:nEnd);

%%% Test for crushes
for k = 1:2:sz
    [d] = p_poly_dist(xMan, yMan, Terrain(k,:), Terrain(k+1,:) )  ;

    if abs(d) < minDistTH
        isCrushed = true;
        objectNum = k;
        break
    end % if
end % for k

%%% Test for new Position if th bat is in the obstacle
if isCrushed
    % [isBatInside,~] = inpolygon(xB, yB, Terrain(objectNum,:), Terrain(objectNum+1,:) )  ;
    isBatInside = any(inpolygon(xMan, yMan, Terrain(objectNum,:), Terrain(objectNum+1,:) )  );
    if objectNum == 1 % the room limits, check if the bat Crossed them and is outside ther polyogon
        RecoverFlag = ~isBatInside;
    else % the bat crushess into one of the Obtacles
        RecoverFlag =   isBatInside;
    end  % if objectNum === 1 %

    %%% id croees the line go back to the Beginning (as Viziny Said ...)
    %%% and turn to the oposite direction
    if RecoverFlag
        xRecover    = BAT.xBati(nStart); % the stating pot of thr Call
        yRecover    = BAT.yBati(nStart);
        [d2, xPoly, yPoly] = p_poly_dist(xRecover, yRecover, Terrain(objectNum,:), Terrain(objectNum+1,:) )  ;
        tetaRecover =  atan2(yRecover- yPoly, xRecover- xPoly);
    else % new 07Jan23
        xRecover   = xB; % the end point of thr Call
        yRecover   = yB;
        [d2, xPoly, yPoly] = p_poly_dist(xRecover, yRecover, Terrain(objectNum,:), Terrain(objectNum+1,:) )  ;
        tetaRecover =  atan2(yRecover- yPoly, xRecover- xPoly);
    end %  if RecoverFlag
end % if isCrushed


debug_flag = 0;
if debug_flag
    figure; hold on
    plot(xB, yB, 'xr')
    plot(xMan, yMan, 'r')
    for k = 1:2:sz
        plot(Terrain(k,:), Terrain(k+1,:),'.k-')
    end % for k
%     if RecoverFlag
        plot(xPoly, yPoly, '*r')
        plot(xRecover, yRecover,'*b')
        plot([xRecover, xRecover+1*cos(tetaRecover)], [yRecover, yRecover+1*sin(tetaRecover)],'^b-', 'LineWidth', 1.5)
%     end
end % if debug_flag

