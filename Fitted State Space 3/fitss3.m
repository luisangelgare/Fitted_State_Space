function fitted_ss = fitss3(Yqdw_full, fd0, poles1, pole_tol, red_mode)
    
    % fitss3: Version 3 of the function that generates a 
    % state-space model using vector fitting

    % Inputs descritption:
    % - Yqdw_full is the [3, 3, fd0] matrix with all the frequency (Hz) measurements.
    % - fd0 is a 1xfd0 real-positive vector which contains the frequency measurements in Hz.
    % - poles1 are the starting complex poles for the fitting process.
    % - pole_tol is the convergence tolerance for the pole iteration.
    % - red_mode is the reduced-order methodology for the fitted state-space:
    % - 1 -> State-space balanced realization (BR) [default].
    % - 2 -> State-space by singular value descomposition (SVC).
    % - 3 -> Full order state-space (without any reduction method)

   %% Validate dimensions of Yqd0_full

    % Validate dimensions of Yqd0_full
    dimensions = size(Yqdw_full);
    if numel(dimensions) ~= 3 || dimensions(1) ~= 3 || dimensions(2) ~= 3 || dimensions(3) ~= length(fd0)
        error('Input Yqd0_full must have dimensions [3, 3, length(fd0)]. Simulation stopped.');
    end


    %% Validate weight

    % - weight is the 1xfd0 matrix that contains the weighted vector (1 by default)
    weight=linspace(1, 1, length(fd0)); % Vector de ponderacion

    %% Set vector fitting options

    opts.relax = 1; % Allow relaxation for non-triviality constraint
    opts.stable = 1; % Force poles stability
    opts.asymp = 2; % Include D and E in the fitting
    opts.skip_pole = 0; % Do not skip pole identification
    opts.skip_res = 0; % Do not skip residue identification
    opts.cmplx_ss = 0; % Use real-valued state-space models
    opts.spy1 = 0; % Disable first fitting stage plot
    opts.spy2 = 0; % Enable magnitude fitting plot
    opts.logx = 1; % Logarithmic x-axis for frequency
    opts.logy = 1; % Logarithmic y-axis for magnitude
    opts.errplot = 1; % Include error deviation in plots
    opts.phaseplot = 0; % Disable phase plot
    opts.legend = 1; % Include legends in plots

     %% Extract submatrices

    % Extract Yqq, Yqd, Yqw, Ydq, Ydd, Ydw, Ywq, Ywd, Yww as 1xfd0 vectors directly
    Yqq = transpose(squeeze(Yqdw_full(1,1,:))); % Maintain dimensions 1xfd0
    Yqd = transpose(squeeze(Yqdw_full(1,2,:))); % Maintain dimensions 1xfd0
    Yqw = transpose(squeeze(Yqdw_full(1,3,:))); % Maintain dimensions 1xfd0
    Ydq = transpose(squeeze(Yqdw_full(2,1,:))); % Maintain dimensions 1xfd0
    Ydd = transpose(squeeze(Yqdw_full(2,2,:))); % Maintain dimensions 1xfd0
    Ydw = transpose(squeeze(Yqdw_full(2,3,:))); % Maintain dimensions 1xfd0
    Ywq = transpose(squeeze(Yqdw_full(3,1,:))); % Maintain dimensions 1xfd0
    Ywd = transpose(squeeze(Yqdw_full(3,2,:))); % Maintain dimensions 1xfd0
    Yww = transpose(squeeze(Yqdw_full(3,3,:))); % Maintain dimensions 1xfd0
    
    %% Perform vector fitting for each subsystem

    % Extract subsystems from the input 2x2xfd0 Yqd0_full array
    % Yqq fitting
    system1 = PCCF1(Yqq, 1i*2*pi*fd0, poles1, weight, pole_tol, red_mode, opts);

    % Yqd fitting
    system2 = PCCF1(Yqd, 1i*2*pi*fd0, poles1, weight, pole_tol, red_mode, opts);

    % Yqw fitting
    system3 = PCCF1(Yqw, 1i*2*pi*fd0, poles1, weight, pole_tol, red_mode, opts);

    % Ydq fitting
    system4 = PCCF1(Ydq, 1i*2*pi*fd0, poles1, weight, pole_tol, red_mode, opts);

    % Ydd fitting
    system5 = PCCF1(Ydd, 1i*2*pi*fd0, poles1, weight, pole_tol, red_mode, opts);

    % Ydw fitting
    system6 = PCCF1(Ydw, 1i*2*pi*fd0, poles1, weight, pole_tol, red_mode, opts);

    % Ywq fitting
    system7 = PCCF1(Ywq, 1i*2*pi*fd0, poles1, weight, pole_tol, red_mode, opts);

    % Ywd fitting
    system8 = PCCF1(Ywd, 1i*2*pi*fd0, poles1, weight, pole_tol, red_mode, opts);

    % Yww fitting
    system9 = PCCF1(Yww, 1i*2*pi*fd0, poles1, weight, pole_tol, red_mode, opts);

    %% Build unified MIMO state-space matrices

    % Construct the A matrix
    Afull = blkdiag(system1.A, system2.A, system3.A, ...
                    system4.A, system5.A, system6.A, ...
                    system7.A, system8.A, system9.A);

    % Construct the B matrix with zero-padding
    Bfull = [system1.B, zeros(size(system1.B)), zeros(size(system1.B));
             zeros(size(system2.B)), system2.B, zeros(size(system2.B));
             zeros(size(system3.B)), zeros(size(system3.B)), system3.B;
             system4.B, zeros(size(system4.B)), zeros(size(system4.B));
             zeros(size(system5.B)), system5.B, zeros(size(system5.B));
             zeros(size(system6.B)), zeros(size(system6.B)), system6.B;
             system7.B, zeros(size(system7.B)), zeros(size(system7.B));
             zeros(size(system8.B)), system8.B, zeros(size(system8.B));
             zeros(size(system9.B)), zeros(size(system9.B)), system9.B];

    % Construct the C matrix with zero-padding
    Cfull = [system1.C, system2.C, system3.C, zeros(size(system4.C)), zeros(size(system5.C)), zeros(size(system6.C)), zeros(size(system7.C)), zeros(size(system8.C)), zeros(size(system9.C));
             zeros(size(system1.C)), zeros(size(system2.C)), zeros(size(system3.C)), system4.C, system5.C, system6.C, zeros(size(system7.C)), zeros(size(system8.C)), zeros(size(system9.C));
             zeros(size(system1.C)), zeros(size(system2.C)), zeros(size(system3.C)), zeros(size(system4.C)), zeros(size(system5.C)), zeros(size(system6.C)),system7.C, system8.C, system9.C];
    
    % Construct the D matrix
    Dfull = [system1.D, system2.D, system3.D;
             system4.D, system5.D, system6.D;
             system7.D, system8.D, system9.D];

    %% Define state-space model

    % No state names provided as per instruction
    fitted_ss = ss(Afull, Bfull, Cfull, Dfull, ...
                   'StateName', {}, 'inputname', {'vq_fit', 'vd_fit','w_ref'}, ...
                   'outputname', {'iq_fit', 'id_fit','w_VSC'});

    fitted_ss = minreal(fitted_ss);

    %% Check system stability
    
    % Stop simulation if the system is unstable
    if ~isstable(fitted_ss)
        error('The fitted state-space model is unstable. Simulation stopped.');
    end
end
