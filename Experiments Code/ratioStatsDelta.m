function [meanZ, stdZ] = ratioStatsDelta(mx, stdx, my, stdy)
% ratioStatsDelta - Approximates mean and std of z = x/y using the delta method
%
% Inputs:
%   mx   - Mean of x
%   stdx - Standard deviation of x
%   my   - Mean of y
%   stdy - Standard deviation of y
%
% Outputs:
%   meanZ - Approximated mean of z = x/y
%   stdZ  - Approximated standard deviation of z = x/y

    % Check that the denominator mean is not zero
    if my == 0
        error('Mean of y (denominator) must not be zero for delta method.');
    end

    % Approximated mean of z
    meanZ = mx / my;

    % Approximated variance of z using delta method
    varZ = (stdx / my)^2 + ((mx * stdy) / my^2)^2;
    stdZ = sqrt(varZ);
end

% % function [meanZ, stdZ] = ratioStats(mx, stdx, my, stdy, N)
% % % ratioStats - Estimates the mean and standard deviation of z = x/y
% % % using Monte Carlo sampling
% % %
% % % Inputs:
% % %   mx   - Mean of x
% % %   stdx - Standard deviation of x
% % %   my   - Mean of y
% % %   stdy - Standard deviation of y
% % %   N    - (Optional) Number of samples [default: 1e6]
% % %
% % % Outputs:
% % %   meanZ - Estimated mean of z = x/y
% % %   stdZ  - Estimated standard deviation of z = x/y
% % 
% %     if nargin < 5
% %         N = 1e6; % Default number of samples
% %     end
% % 
% %     % Generate samples
% %     x = mx + stdx .* randn(N, 1);
% %     y = my + stdy .* randn(N, 1);
% % 
% %     % Remove problematic divisions
% %     valid = abs(y) > 1e-6;
% %     z = x(valid) ./ y(valid);
% % 
% %     % Compute statistics
% %     meanZ = mean(z);
% %     stdZ = std(z);
% % end
