function fitted_ss = MIMOn4sid(Yqd_full, fd0, Ts, nx)
    % fit_mimo_n4sid - Fit a unified MIMO state-space model using N4SID algorithm
    %
    % Inputs:
    %   Yqd_full : 3x3xN matrix containing frequency response data
    %   fd0      : Vector of frequency points (in Hz)
    %   Ts       : Sampling time (0 for continuous systems)
    %   nx       : Desired model order
    %
    % Output:
    %   fitted_ss : Unified state-space model (ss object)

    % Check input dimensions
    if ndims(Yqd_full) ~= 3
        error('Input Yqd_full must be a 3x3xN 3D matrix.');
    end
    if size(Yqd_full,1) ~= 3 || size(Yqd_full,2) ~= 3
        error('Input Yqd_full must have dimensions 3x3xN.');
    end
    if size(Yqd_full,3) ~= length(fd0)
        error('Third dimension of Yqd_full must match the length of fd0.');
    end

    % Initialize cell array to store individual SISO models
    systems = cell(3,3);

    % Fit each SISO model using n4sid
    for i = 1:3
        for j = 1:3
            resp = squeeze(Yqd_full(i,j,:));
            Gij = idfrd(resp, fd0, Ts, 'Units', 'Hz');
            systems{i,j} = n4sid(Gij, nx);
        end
    end

    % Extract individual models
    system1 = systems{1,1}; system2 = systems{1,2}; system3 = systems{1,3};
    system4 = systems{2,1}; system5 = systems{2,2}; system6 = systems{2,3};
    system7 = systems{3,1}; system8 = systems{3,2}; system9 = systems{3,3};

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
             zeros(size(system1.C)), zeros(size(system2.C)), zeros(size(system3.C)), zeros(size(system4.C)), zeros(size(system5.C)), zeros(size(system6.C)), system7.C, system8.C, system9.C];

    % Construct the D matrix
    Dfull = [system1.D, system2.D, system3.D;
             system4.D, system5.D, system6.D;
             system7.D, system8.D, system9.D];

    % Create the unified state-space model
    fitted_ss = ss(Afull, Bfull, Cfull, Dfull, Ts, ...
                   'inputname', {'vq_fit', 'vd_fit','w_ref'}, ...
                   'outputname', {'iq_fit', 'id_fit','w_VSC'});

    % Check stability of the model
    if ~isstable(fitted_ss)
        warning('The fitted state-space model is not stable.');
    end
end
