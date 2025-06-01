function [Directivity] = BeamDirectivity(Freq, Tetas, TransOrRecFlag)

% calucte the dircetivty of the beam for a given anlge
% using the pattern :
% D(teta) = (2 · J1 (k ª a ª sin(?)) / k ª a ª sin(?))^2
% Inputs:
%   Freq - the frequency of the signal
%   Tetas- the angles to calculte
%   TransOrRecFlag =Transmit or Recieved - change the gain by the area of
%   the surface ("mouth" /  "ear")
%
% Correct Bug Sep2022 - Abs : Alpha = abs(K0.*AntRadius.*sin(Tetas));

MouthR = 0.003;% 0.003; 0.006
EarR= 0.007; %0.007; 0.012

switch TransOrRecFlag
    case 'Transmit'
        %%%% For Transmit Gain =1 ( 0 dB)
        AntRadius = MouthR; % [m] - the mouth of the bat
        Gain = 1;
        
    case 'Recieved'
        AntRadius = EarR; % [m] - the mouth of the bat
        AntArea = AntRadius^2*pi;
        mSoundVelocity = 343; %m/sec
        Lambda = mSoundVelocity ./ (Freq*1e3); % m
        K0 = 2*pi ./ Lambda;
        Gain = 4*pi.*AntArea ./ Lambda.^2;
end %switch TransOrRecFlag

%%%% For Transmit Gain =1 ( 0 dB)
% AntArea = AntRadius^2*pi;
mSoundVelocity = 343; %m/sec
Lambda = mSoundVelocity ./ (Freq*1e3); % m
K0 = 2*pi ./ Lambda;
% Gain = 4*pi.*AntArea ./ Lambda.^2;


Alpha = abs(K0.*AntRadius.*sin(Tetas));
% Dealing with divide by zero
Alpha(Alpha==0) = 1e-10;
% 2Pi directivity
Directivity = Gain .* abs( ((2*besselj(1,Alpha) ./ Alpha) ) ).^2; %%% adding sqr 28/8/8/18
% Directivity = Gain .* abs( ((2*sin(Alpha) ./ Alpha) ) );

%%% Adding 0-20dB att to back angles
B = 1.99; A = 1.98/pi;
Mask = ones(size(Tetas));
Ind = find(abs(Tetas)> pi/2);
Mask(Ind) = -A*abs(Tetas(Ind)) + B;
Directivity = Directivity .* Mask;

%Dircetivity = Gain .* abs( ((2*besselj(1,K0.*AntRadius.*sin(Tetas)) ./ Alpha)) );

