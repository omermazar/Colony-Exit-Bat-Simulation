function [classifierTestVec, results] = PreyClassifyFilterBank(sig_echo, FBResponseIn, nTimesToCheck, sig_target, FBResponseTargetRef, tstart , fs, fc, plotFlag)

%% Classify detected peaks - whether prey or not
% classifierTestVec = PreyClassifyFilterBank(sig_echo, FBResponseIn, detectedPksToCheck, sig_target, FBResponseTargetRef, tstart , fs, fc)
% Returns a logical vector of with the classifier answer ('0'/'1') for each of tested peaks defined detectedPksToCheck. 
% 
% [classifierTestVec, classifierResults] = PreyClassifyFilterBank(___) 
% also returns classifierResults, that is a struct with the full results of
% the calssifiier tests
%
% INPUTS:
%   sig_echo - the signal to be tested (samples), padded with zeros in the
%   start and end of the vector
%   FBResponseIn - the FilterBank Response to sig_echo
%   nTimesToCheck - the times (samples) of the required segment to check
%           in the Inpout Vector
%   sig_target - the reference signature of the target
%   FBResponseTargetRef - the FilterBank response to sig_target
%   tstart - the time of the starting and anding paddings
%   fs - the sample frequency of the signals
%   fc - the center frequencies of the FilterBank

%%%% 

%% Parameters and Init
% tstart = 1/fs;
response_time = 0.6e-3;
nresponse = round(response_time * fs);
% plotFlag  = 1; 
nstart    = tstart * fs;
nPeaks    = numel(nTimesToCheck);
classifierTestVec = false(1, nPeaks);
classifierResults = cell(nPeaks,1);
results = struct();

% methods of the classifier
methods.target_flg        = 0;  % for refence parmeters define as 1
methods.normalzation      = 1;  % noramlize the signals according  to their energy
methods.standarzization   = 0;  % standardization of the singnals (mean =1; std =0), (Doesnt work well) Should be standarzization or  normalzation
methods.correaltionTest   = 1;  % test correlation coefficent between the signals%
methods.pksTrghsTest      = 0;  % test correlation coefficent between the signals%
methods.mseTest           = 0;  % test correlation coefficent between the signals%
methods.cochleogramMatrix = 1;  % run calssifier on the time-frequency response matrix
methods.frequencyResponse = 1;  % run calssifier on the frequency response vector
methods.timeResponse      = 0;  % run calssifier on the time response vector

%% The reference Signature
%%%% Active Channels
active_ind = FBResponseTargetRef.Active_Channels;
n_active   = sum(active_ind);
f = fc(active_ind);

%%%% reference Matrix and Frequency Response (FR)
mat_target_full =  reshape([FBResponseTargetRef.Channels(active_ind).shifted_response], [], n_active);
ix_mat_target   = round((nstart-nresponse):(nstart+nresponse));
mat_target      = mat_target_full( ix_mat_target, : );
ix_start_target = round(nresponse+1);
FR_target       = mat_target(ix_start_target, : );

% [~, pks_target, time_pks_target] = analyze_echo2ref(sig_target, sig_target, mat_target, FR_target, mat_target, FR_target, true, [], [], f, fs);
if methods.pksTrghsTest
    methods.target_flg  = 1;
    [~, pks_target]  = ClassifierCompareToRef( FR_target, FR_target, methods);
else
    pks_target = [];
end % if methods.pksTrghsTest
methods.target_flg  = 0;

 if plotFlag
    fig = figure;
    
    fig.Position = [ 200 100 560 600];
    ax1 = subplot(3,2,1);  % title(ax1, 'Target')
    ax2 = subplot(3,2,2); % title(ax2, 'detection' )
    ax3 = subplot(3,2,3); hold on;  % title(ax3, 'Frequncy Response, target (Black)')
    ax4 = subplot(3,2,5); hold on; grid minor % title(ax3, The signal')
    
    ax3.Position(3) = ax2.Position(1) + ax2.Position(3) - ax1.Position(1);
    ax4.Position(3) = ax2.Position(1) + ax2.Position(3) - ax1.Position(1);

    imagesc(ax1, mat_target');
    title(ax1, 'Target')
    
    plot(ax3, f, FR_target/norm(FR_target), 'k.-', 'LineWidth', 1.5 ,'DisplayName', 'Target Signalture')
    title(ax3, 'Frequncy Response, target (Black)') 
    
    plot(ax4, 10*log10(abs(sig_echo)) )
    plot(ax4, ax4.XLim, zeros(1,2) )
    ax4.YLim = [-20, max(10*log10(abs(sig_echo))+10)];
    
    title(ax4, 'The Full Echo')
    sgtitle('Classifier')
 end % if plotFLag

%% The Recieved Echo to Check
mat_full = reshape([FBResponseIn.Channels(active_ind).shifted_response], [], n_active);
% pad the start and end periods with zeros
mat_full = padarray(mat_full, nstart);
sig_echo = padarray(sig_echo', nstart)';

%%%% check the claiffier in each detection of the Detector of the input signa 
detectednTimes = nstart + nTimesToCheck; % FBResponseIn.det_ipks(nTimesToCheck); %  fix the detection times according to the added zeros 
k = 1;
for inTime = detectednTimes
    % cut the signal in the detection times
    ix_curDetection = round((inTime-nresponse):(inTime+nresponse));
    mat_currEcho    = mat_full(ix_curDetection, : );

    % the Frequency respons is taken from the detection time
    FR_currEcho  = mat_currEcho(round(nresponse+1), : ); 
    max_ix       = min(numel(sig_echo), inTime + numel(sig_target) + nresponse );
    sig_currEcho = sig_echo( (inTime-nresponse):max_ix );
    
    % the result
%     [classifierResultsXXX, pks_out, time_pks_out, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
%         mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2]  = ...
%         analyze_echo2ref(sig_target, sig_currEcho, mat_target, FR_target, mat_currEcho, FR_currEcho, ...
%         false, pks_target, time_pks_target, f, fs );

    classifierResults{k}       = ClassifierCompareToRef( FR_target, FR_currEcho, methods, pks_target, mat_currEcho, mat_target);
    classifierResults{k}.ntime = nTimesToCheck(k);

    classifierTestVec(k) = ClassifierDecision(classifierResults{k}, methods);
    classifierResults{k}.isClassified = classifierTestVec(k);

    if plotFlag
        imagesc(ax2, mat_currEcho');
        title(ax2, ['Echo Detection: ', num2str(k)])
        currIx = (inTime-nresponse):(max_ix);

        plot(ax4, currIx, sig_currEcho)
        if classifierTestVec(k)
            txtLegend = ['Echo:', num2str(k), '+'];
            plot(ax4, inTime, ax4.YLim(2)-10, '+', 'MarkerSize', 6, 'LineWidth', 2)
        else
            txtLegend = ['Echo:', num2str(k)];
            plot(ax4, inTime, ax4.YLim(2)-10, 'o' ,'MarkerSize', 6, 'LineWidth', 2)
        end 
        plot(ax3, f, FR_currEcho/ norm(FR_currEcho), '.-', 'DisplayName', txtLegend)
        hold(ax3, 'on')
        legend(ax3,'Location','eastoutside')
        plot(ax3, f, FR_target/norm(FR_target), 'k.-', 'LineWidth', 1.5 ,'DisplayName', 'Target Signalture')
%         hold(ax3, 'off')
%        pause;
        
    end % if plotFLag
    k = k+1;
end % for inTime = detectednTimes

%% Output
if ~isempty(nTimesToCheck)
    fn = fieldnames(classifierResults{1});
    tempVec = nan(1,numel(classifierResults));

    for n = 1:numel(fn)
        for k =1:numel(classifierResults)
            tempVec(k) = classifierResults{k}.(fn{n});
        end % f k
        results.(fn{n}) = tempVec;
    end % for n
end % if ~isempty(nTimesToCheck)



% if plotFlag
%     myClassifierResultsPlot(results)
% end





