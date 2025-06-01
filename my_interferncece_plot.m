function [Masking_strct, fig] = my_interferncece_plot(BatDATA, BatNum, percentile) 

plot_flg = 0;
if nargin < 3
    prc_txt = repmat(' ',size(BatNum'));
else
    prc_txt = string(num2str(percentile'));
end % if nargin
if nargin ==1
    BatNum =1;
end % if nargin

n_plots = numel(BatNum);
ts = BatDATA.AllParams.SimParams.SampleTime;
MaxTime = BatDATA.AllParams.SimParams.SimulationTime;
tt = [0:ts:MaxTime]; % Time Vector

Masking_strct.percentile = percentile;
Masking_strct.Bats = BatNum;
Masking_strct.time_vec = tt;

if BatDATA.AllParams.SimParams.swarm_initial_teta_std == 0
    rand_txt = '';
else
    rand_txt = ' random angle';
end

if plot_flg
    fig = figure;
end % plot_flg

for k=1:n_plots
    
    Masking_strct.Bat(k).BatNum = BatNum(k);
    Masking_strct.Bat(k).Percentile= percentile(k);
    Masking_strct.Bat(k).own_calls_lvl = BatDATA.BAT(BatNum(k)).BatSonarEchosMat(1,:);
    Masking_strct.Bat(k).masking_lvl = BatDATA.BAT(BatNum(k)).AllInterPulses;
    Masking_strct.Bat(k).rx_lvl_from_conspecific = BatDATA.BAT(BatNum(k)).Consps_EchosVec;
    
    if plot_flg
        ax = subplot(n_plots, 1, k);
        hold on
        plot(ax,tt, 10*log10(abs(BatDATA.BAT(BatNum(k)).AllInterPulses)), 'r-', 'linewidth',1, 'DisplayName','Masking' )
        plot(ax, tt, 10*log10(abs(BatDATA.BAT(BatNum(k)).BatSonarEchosMat(1,:) )) , 'k','linewidth',0.5,  'DisplayName','Own Calls' )
        plot(ax, tt, 10*log10(abs(BatDATA.BAT(BatNum(k)).Consps_EchosVec)), 'g','linewidth',1, 'DisplayName','Cons. Echoes' )
        ax.YTick = [0:20:110];
        grid on
        
        ylabel('level(dB-SPL)')
        posx = num2str(BatDATA.BAT(BatNum(k)).xBatPos(1));
        posy = num2str(BatDATA.BAT(BatNum(k)).yBatPos(1));
        title(['Percentile: ', char(prc_txt(k)) ,'%', ', StartPos:', posx,',',posy ]);
    end % plot_flg
end % for k

if plot_flg
    legend
    xlabel('time(sec)')
    
    suptitle(['Total Bats: ', num2str(BatDATA.AllParams.SimParams.TotalBatsNumber), ...
        ', std: ', num2str(BatDATA.AllParams.SimParams.swarm_initial_ybat_std), 'm', rand_txt])
    
    set(fig,'Units','centimeters','Position',[ 6 3 21 18])
end % plot_flg