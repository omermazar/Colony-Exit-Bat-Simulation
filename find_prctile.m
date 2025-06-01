function [Masking_strct, bat_nums] = find_prctile(BatDATA, percentile)

% returns the batc in the requried percentiles
if nargin < 2
    percentile = [50, 80, 90];
end % nargin

nsamples= size(BatDATA.BAT(1).xBatPos,2);
bat_nums = zeros(size(percentile));
all_y = [BatDATA.BAT.yBatPos];
all_y = all_y(1:nsamples:end);
for k=1:numel(percentile)
    py = prctile(all_y, percentile(k));
    [ d, ix ] = min( abs( all_y- py ) );
    bat_nums(k) = ix;
end

Masking_strct = my_interferncece_plot(BatDATA, bat_nums, percentile) ;
