%% test rx level Script

%% Gain
  

%% Attenuation
Tetas = 0;
Distances = 0.05:0.1:10; %1; % 0.05:0.1:5;
TargetArea = 0.1; % 0.15^2*pi; %0.02^2*pi; % 0.0025; % 0.02^2*pi;
Freq = 40;% 25; %[20:10:70]; % [1,2,5,10:10:70];
% Sonar
[sonarAttenuation] = CalculateSignalAtten(Distances, Tetas, TargetArea, Freq, 'meters');

% Target Strength
SigmaRCS = 1;
Lambda = 343./(Freq*1e3);
rTarget =  sqrt(TargetArea/pi);
TargetRCS = TargetArea*ones(size(Lambda));
ixSmall = rTarget./Lambda < 1.5;
%     ixBig = rTarget./Lambda > 4;
TargetRCS(ixSmall) = SigmaRCS*4*pi*TargetArea^2 ./ Lambda(ixSmall).^2; % New Nov22
%     TargetRCS(~ixSmall & ~ixBig) = SigmaRCS*2*pi*TargetArea; % New Nov22
%     TargetRCS(ixBig) = TargetArea; % New Nov22
TS = 10*log10(TargetRCS/(4*pi))

sphereRCS = TargetArea; % New Nov22
planarRCS = 4*pi*TargetArea^2 ./ Lambda.^2; % New Nov22
maxSphere =  4*pi*TargetArea;
figure;
hold on
plot(Freq, 10*log10(TargetRCS/(4*pi)), '*-', 'DisplayName', 'Full')
plot(Freq, 10*log10(sphereRCS/(4*pi)).*ones(size(Freq)), '+-', 'DisplayName', 'sphere')
plot(Freq, 10*log10(planarRCS/(4*pi)), 'o-')
plot(Freq, 10*log10(maxSphere/(4*pi).*ones(size(Freq))), '--')
% TargetRCS = TargetArea;

% mask
BatMouthGain = BeamDirectivity(Freq, Tetas, 'Transmit');
BatEarGain = BeamDirectivity(Freq, Tetas, 'Recieved');
AlphaAtmAtt = 0.038*Freq - 0.3; 
Lambda = 343/(Freq*1e3);
maskAttenuation =  BatMouthGain .* BatEarGain .* ...
    10.^(-1.*AlphaAtmAtt./10.*(Distances-0.1)) .* (Lambda.^2  ./ (4*pi.*Distances).^2 );

% wrong calculation
%%% change ReflectorType = 'Planar'
[sonarAttenuationXXX] = CalculateSignalAtten(Distances, Tetas, TargetArea, Freq, 'meters');
% Distances = 0.05:0.1:5;
% Tetas = 0;
figure
hold on
plot(Distances, 10*log10(abs(sonarAttenuation)), 'DisplayName', 'Sonar Echo')
plot(Distances, 10*log10(abs(maskAttenuation)), 'DisplayName', 'maskConsp')
xlabel('Disrtance (m)')
ylabel('Attenuation (dB)')
grid
legend
plot(Distances, 10*log10(abs(sonarAttenuationXXX)), 'k-', 'DisplayName', 'Wrong Sonar')


%% Target Area
r = [0.005:0.005:0.2]; % the radious of the target
PlanarRCS = 4*pi*(pi*r.^2).^2 ./ Lambda.^2;
SphereRCS = (pi*r.^2).^2 ;
figure
hold on
plot(r, PlanarRCS, 'DisplayName', 'Planar')
plot(r, SphereRCS, 'DisplayName', 'Sphere')
legend
xlabel('Target radius (m)')
ylabel('RCS (m^2)')

% Target strength (TS) is equal to 10 log10(σbs/(1 m2)) dB, where σbs is the differential backscattering cross section. Backscattering cross section is 4πσbs.