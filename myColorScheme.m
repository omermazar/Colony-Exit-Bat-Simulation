function [colStruct, colObj] = myColorScheme(Mode, nBats, nPreys, Object, plotFlag)

%% Input: Mode 'foraging' (default) | 'swarm' | 'caveExit'
%        Object
if nargin < 5
    plotFlag = false;
end

if nargin < 4
    Object = 'Prey';
end
 
if nargin < 3
    nPreys = 10;
end

if nargin < 2
    nBats = 10;
end


if nargin < 1
    Mode = 'foraging';
end

%% The colormap
%%%% Acoustics
colStruct.ownCalls           = [0 0 0]; % black
colStruct.Prey               = [0.4660    0.6740    0.1880]; % green
colStruct.conspsEchoes       = [0    0.4470    0.7410]; % blue
colStruct.clutter            = [0.9290    0.6940    0.1250]; % yellow
colStruct.conspsCalls        = [0.9, 0.3, 0]; % [1 0 0]; %[0.8500    0.3250    0.0980] ; % red
colStruct.conspsObsEchoes    = [1 0.4 0]; % 0.7*[ 1 1 1 ]; % grey
colStruct.conspsConspsEchoes = [0.8 0 0.2]; % [0.5928    0.2208    0.6672]; % purple
colStruct.MaskingSignals     = [1 0 0]; % red
colStruct.wantedSignals      = [0 1 0]; % green
%%%% XY and Animations
% objects
colStruct.Terrain        = [0 , 0 , 1];
colStruct.BatPos         = 0.4*ones(nBats,3); %grey % lines(nBats);
colStruct.SelectedBatPos = [0.9    0.9    0.250]; % yelow
colStruct.PreysPos       = lines(nPreys);
colStruct.EcholocatePos  = [0.3010    0.7450    0.9330]; % cyan
% manuvers
colStruct.ForagingMan = [0 0 0]; %Black
colStruct.Approach    = [ 0 1 0];  % Green
colStruct.BuzzMan     = [ 0.2 0.2 0.2 ]; % Grey
colStruct.ObsMan      = [0.7 0.7 0]; % yelow dark
colStruct.CaveExit    = [0 1 1]; %% cyan
colStruct.AvoidBatMan = [1 1 0]; % Yellow;
colStruct.FollowWall  = [ 0 1 0];  % Green sane as Approach

% detections
colStruct.preyFinding   = [0 1 0]; % green
colStruct.conspFinding  = [0 0.8 0]; % green
colStruct.obsFinding    = [0 0.5 0]; % green
colStruct.preyEstimate  = [0 1 0]; % green
colStruct.conspEstimate = [0 0.8 0]; % green
colStruct.obsEstimate   = [0 0.5 0]; % green

%%%  events
colStruct.catchPrey = [ 0 0 0]; % black
colStruct.crushes   = [0.6 0.6 0.6] ; % grey
colStruct.detections= [0 1 0] ; % gren
% maskig color
colStruct.masking = [1 0 0]; %red

switch Mode
    case 'foraging' 
       colStruct.Prey         = [0.4660    0.6740    0.1880]; % green
       colStruct.conspsEchoes = [0    0.4470    0.7410]; % blue
       colStruct.clutter      = [0.9290    0.6940    0.1250]; % yellow
        
    case'swarm' 
       colStruct.Prey         = [0    0.4470    0.7410]; % blue
       colStruct.conspsEchoes = [0.4660    0.6740    0.1880]; % green
       colStruct.clutter      = [0.9290    0.6940    0.1250]; % yellow

    case 'caveExit'
       colStruct.Prey         = [0.9290    0.6940    0.1250]; % yellow
       colStruct.conspsEchoes = [0    0.4470    0.7410]; % blue
       colStruct.clutter      = [0.4660    0.6740    0.1880]; % green
end % switch

if nargout >= 2
    colObj = colStruct.(Object);
end

%% plot legend
if plotFlag || nargout == 0
    allCols = fields(colStruct);
    fig = figure;
    hold on
    nCols = numel(allCols);
    n=0;
    for k = 1: nCols
        if ~ismember(allCols{k}, {'BatPos', 'PreysPos'} )
            n=n+1;
            plot([1, 1.005], ones(1,2)*(nCols-n+1), 'o', 'MarkerSize', 9,  'MarkerFaceColor', colStruct.(allCols{k}), 'MarkerEdgeColor', 'none' )
            text(1.05, 1*(nCols-n+1), allCols{k}, 'FontSize', 8)
            xlim([0.98, 1.3])
            fig.Position = [88.0000  242.0000  200  500];
            title('BatSimulation Colors')
        end %
    end % for k
end % if plot flaf