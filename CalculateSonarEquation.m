function[Attenuation] = CalculateSonarEquation(Distance, Angle, TargetRCS, Freq, ...
        BatFindEnvironment(xCurrent,yCurrent,AngleView,BeamWidth, ...
        DetectionDistance, Terrain, PreysPositions)
    
    
    %%% CONSTS
    AttenuationConst = 1; % Just parameter
%     AlphaAtmAtt = 0.23; %   Atmosphere attenuation[m^(-1)] ...
    %  1dB loos/meter % 0.3 for 40kHa (at 30C,70% humidty],alpha =  0.1 for 28kH
    AlphaAtmAtt = 0.038*Freq - 0.3; % dB/m from table, @20Cel.50% humidity Atmosphere attenuation[m^(-1)] ...
                    %  1dB loose/meter % 0.3 for 40kHa (at 30C,70% humidty],alpha =  0.1 for 28kH
                    
    BatMouthGain = BeamDirectivity(Freq, Teta, 'Transmit');
    BatEarGain = BeamDirectivity(Freq, Teta, 'Recieved');
    
    Attenuation = AttenuationConst * BatMouthGain * BatEarGain * TargetRCS * ...
        exp(-2*AlphaAtmAtt.*(xyDistances-0.1)) ./ (4*pi*xyDistances).^4 .*...
        exp(1i*4*pi.*xyDistances./lambda);