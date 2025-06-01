function[Attenuation] = CalculateSignalAtten(Distances, Tetas, TargetArea, Freq, CoordinationFlag, varargin)
%  
% function[Attenuation] = CalculateSignalAtten(Distances, Angles,
% TargetRCS, Freq, CoordinationFlag,  (SoundVelocirty,, RefelectorTYpe))
% calculate te attenuatian of the reciveid Sonar Signal
% Inputs - 
%       Distances- vector of the distances from target
%       Tetas- the angle between the bat and the prey
%       Target RCS - the target cross section
%       Freq- the frequency [kHz]
%       CoordinationFlag = 'meters' or 'xy'
%       if CoordinationFlag = 'xy' than xyresolution is neede in varagin{1}
%       varagin{2} (optional) - AllParams
%       varagin{3} (optional) - 
% Outputs -
%       Attenuation - vector of the attenuations

%%%% Inputs
switch CoordinationFlag
    case 'meters'
        mDistances = Distances+0.005; % avoid deviding by zero
    case 'xy'
        xyResolution= varargin{1};
        mDistances = (Distances+0.5)*xyResolution;
end %switch CoordinationFlag

% InputsNum = nargin;
if nargin >= 6
   mSoundVelocity = varargin{1};
else
    mSoundVelocity = 343;  %m/sec
end % if narargin

DirectivityMode = 1;
if nargin >= 7 % DirectivityMode
    AllParams = varargin{2};
    if  any(strcmp('DirectivityMode',fieldnames(varargin{2}.BatSonarParams)))
        DirectivityMode = AllParams.BatSonarParams.DirectivityMode;
    end % if  any(strcmp('DirectivityMode',fieldnames(varargin{2}.BatSonarParams)))
else
   DirectivityMode = 1;  % calculate direcetivi
end % if narargin >=7

if nargin >= 8
   ObsOrPrey = varargin{3};
else
    ObsOrPrey = 'Regular';  %m/sec
end % if narargin

if nargin >= 9
    
   ReflectorType = varargin{4};
else
    ReflectorType = 'Planar'; % 'Sphere'; %'Planar';  %m/sec
end % if narargin

switch ObsOrPrey
    case 'Obs'
        AttenuationConst = 10.^(AllParams.BatSonarParams.ClutterGainDB/10); % Just parameter
    otherwise
        AttenuationConst = 1; % Just parameter
end
% AttenuationConst = 1; % Just parameter

AlphaAtmAtt = 0.038*Freq - 0.3; % dB/m from table, @20Cel.50% humidity Atmosphere attenuation[m^(-1)] ...
                    %  1dB loose/meter % 0.3 for 40kHa (at 30C,70% humidty],alpha =  0.1 for 28kH
Lambda = mSoundVelocity ./ (Freq*1e3); % m
 
switch ReflectorType
    case 'Sphere'
        SigmaRCS = 1;
        TargetRCS = SigmaRCS*TargetArea;
    case 'Planar'
        SigmaRCS = 1;
        rTarget =  sqrt(TargetArea/pi);
        TargetRCS = TargetArea*ones(size(Lambda));
        ixSmall = rTarget./Lambda < 1.5;
        %     ixBig = rTarget./Lambda > 4;
        TargetRCS(ixSmall) = SigmaRCS*4*pi*TargetArea^2 ./ Lambda(ixSmall).^2; % New Nov22
        %     TargetRCS(~ixSmall & ~ixBig) = SigmaRCS*2*pi*TargetArea; % New Nov22
        %     TargetRCS(ixBig) = TargetArea; % New Nov22

end % ReflectorType

BatMouthGain = BeamDirectivity(Freq, Tetas, 'Transmit');

if DirectivityMode
%     BatMouthGain = BeamDirectivity(Freq, Tetas, 'Transmit');
    BatEarGain = BeamDirectivity(Freq, Tetas, 'Recieved');
else %if DirectivityMode
%     BatMouthGain = 1;
    BatEarGain = 1;
end % if DirectivityMode

%%% RADAR equation%%%
Attenuation = AttenuationConst .* BatMouthGain .* BatEarGain .* Lambda.^2 .* TargetRCS .* ...
    10.^(-2*AlphaAtmAtt./10.*(mDistances-0.1)) ./ ( (4*pi).^3.*mDistances.^4 ); % .*...
       % exp(1i*4*pi.*mDistances./Lambda);
    
    