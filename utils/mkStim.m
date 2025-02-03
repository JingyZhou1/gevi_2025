% make temporal stimuli

function stim = mkStim(condition, nframes)

%condition = 'TF';

if ~exist('nframes', 'var'), nframes = 140; end

switch condition
    case 'TF' % temporal frequency
        % place holder for all conditions
        stim = zeros(nframes, 16);
        start_frame = 13; % 10 frames before stimulus starts, I made 2 extra frames for stim. to reach the brain

        % first two conditions are blanks

        % the third condition, starting from frame 11, 20 frames on, 30 off
        % and repeated twice
        stim = mk_cond(stim, start_frame, 20, 30, 2, [3, 10]);

        % the fourth condition, starting from frame 11, 10 frames on, 15
        % frames off, repeated 4 times
        stim = mk_cond(stim, start_frame, 10, 15, 4, [4, 11]);

        % the fifth condition, starting from frame 11, 4 frames on, and 8
        % frames off, repeated 8 times
        stim = mk_cond(stim, start_frame, 4, 8, 8, [5, 12]);

        % the sixth condition, starting from frame 11, 4 frames on, 6
        % frames off, repeated 10 times
        stim = mk_cond(stim, start_frame, 4, 6, 10, [6, 13]);

        % the seventh condition, starting from frame 11, 2 frames on, 4
        % frames off, repeated 16 times
        stim = mk_cond(stim, start_frame, 2, 4, 16, [7, 14]);

        % the eighth condition, 2 frames on, 3 frames off, repeated 20
        % times
        stim = mk_cond(stim, start_frame, 2, 3, 20, [8, 15]);

        % the nineth condition, 1 frame on, 2 frames off, repeated 33 times
        stim = mk_cond(stim, start_frame, 1, 2, 33, [9, 16]);

    case 'frequency_plot'
        stim = [];
        count = 1;

        for k = 1 :  50

            frequency = k; % Set the frequency of the square wave (in Hz)
            sampling_frequency = 100; % Set the sampling frequency (in Hz)
            duration = 3; % Set the duration of the square wave (in seconds)

            % Generate time vector
            t = 0:1/sampling_frequency:duration - 0.01;

            % Generate square wave
            stim(count, :) = square(2 * pi * frequency * t);
            count = count + 1;
         
        end
        stim = reshape(stim', [], 1);
        stim(stim<= 0) = 0;

    case 'contrast'
        ctrst = [1, 0.5, 0.25, 0.125, 0.0625, 0.03125];
        start_frame = 13; % 10 frames before stimulus starts, I made 2 extra frames for stim. to reach the brain
        on_period = 20;
        off_period = 30;
        repeat = 2;

        stim  = zeros(nframes, 14);

        for k1 = 3 : 8
            stim = mk_cond(stim, start_frame, on_period, off_period, repeat, [k1, k1 + 6]);
            stim(:, [k1, k1 + 6]) = stim(:, [k1, k1 + 6]) * ctrst(k1 - 2);
        end

    case 'contrast1'
        ctrst = [1, 0.5, 0.25, 0.125, 0.0625, 0.03125];%[0.03, 0.04, 0.06, 0.125, 0.25, 0.5, 1];
        start_frame = 13;  % 10 frames before stimulus starts, I made 2 extra frames for stim. to reach the brain
        on_period = 6; %20;
        off_period = 14; %30;
        repeat = 5; %3;

        stim  = zeros(132, 6); % for 132 frame, originally

        for k1 = 3 : 8
            stim = mk_cond(stim, start_frame, on_period, off_period, repeat, k1);
            stim(:, k1) = stim(:, k1) * ctrst(k1 - 2);
        end

    otherwise
        disp('Incorrect condition name.')

end



% sub-function that make each temporal frequency condition --------
    function stim = mk_cond(stim, start_frame, n_on, n_off, n_repeats, which_cond)
        count = start_frame;
        for k = 1 : n_repeats
            stim([count: count + n_on], which_cond) = 1;
            count = count + n_on + n_off;
        end

    end



end