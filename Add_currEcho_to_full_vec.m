function [prev_full] = Add_currEcho_to_full_vec(prev_full, curr_vec, start_idx, nSamples)

% the indices for the assignment
Acoustic_idx = start_idx : start_idx + nSamples-1;
AcousticNumel = numel(prev_full);
% full_vec = prev_full;

% check for the end of simulation
if start_idx < AcousticNumel
    
    if Acoustic_idx(end) > AcousticNumel
        Acoustic_idx = Acoustic_idx(1):AcousticNumel;
        curr_vec = curr_vec(1:numel(Acoustic_idx));
    end % if Acoutic_idx(end) > AcousticNumel
    %%% add the result to the acoustic vector
%     full_vec(Acoustic_idx) = prev_full(Acoustic_idx) + curr_vec;
        prev_full(Acoustic_idx) = prev_full(Acoustic_idx) + curr_vec;

else % if start_idx < AcousticNumel
    % do nothing
end % if start_idx < AcousticNumel

                        

