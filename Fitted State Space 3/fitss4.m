function fitted_ss = fitss4(Yqdw_full, fd0, freq_band, pole_tol, red_mode)

    %% Validación de dimensiones
    dimensions = size(Yqdw_full);
    if numel(dimensions) ~= 3 || dimensions(1) ~= 3 || dimensions(2) ~= 3 || dimensions(3) ~= length(fd0)
        error('Input Yqd0_full must have dimensions [3, 3, length(fd0)]. Simulation stopped.');
    end

    %% Vector de ponderación
    weight = linspace(1, 1, length(fd0));

    %% Opciones de vector fitting
    opts.relax = 1; opts.stable = 1; opts.asymp = 2; opts.skip_pole = 0;
    opts.skip_res = 0; opts.cmplx_ss = 0; opts.spy1 = 0; opts.spy2 = 0;
    opts.logx = 1; opts.logy = 1; opts.errplot = 1; opts.phaseplot = 0;
    opts.legend = 1;

    %% Extraer submatrices
    Yqq = transpose(squeeze(Yqdw_full(1,1,:)));
    Yqd = transpose(squeeze(Yqdw_full(1,2,:)));
    Yqw = transpose(squeeze(Yqdw_full(1,3,:)));
    Ydq = transpose(squeeze(Yqdw_full(2,1,:)));
    Ydd = transpose(squeeze(Yqdw_full(2,2,:)));
    Ydw = transpose(squeeze(Yqdw_full(2,3,:)));
    Ywq = transpose(squeeze(Yqdw_full(3,1,:)));
    Ywd = transpose(squeeze(Yqdw_full(3,2,:)));
    Yww = transpose(squeeze(Yqdw_full(3,3,:)));

    %% Generar polos iniciales por subsistema con iteración fija
    
    jw = 1i * 2 * pi * fd0;
    N_init = 2; % Número inicial de polos por subsistema
    subsystems = {Yqq, Yqd, Yqw, Ydq, Ydd, Ydw, Ywq, Ywd, Yww};
    all_poles = [];
    
    for k = 1:9
        % Inicializar polos para el subsistema k
        poles_k = -100 - 1i * linspace(2*pi*freq_band(1), 2*pi*freq_band(2), N_init);
    
        % Iterar vectfit3 al menos 8 veces para mejorar precisión
        for iter = 1:18
            [~, new_poles, ~, ~] = vectfit3(subsystems{k}, jw, poles_k, weight, opts);
            poles_k = new_poles; % Actualizar polos para la siguiente iteración
        end
    
        % Almacenar los polos finales del subsistema k
        all_poles = [all_poles; poles_k(:)];
    end


    %% Colapsar polos repetidos
    collapsed_poles = [];
    tolerance = 0.03; % Tolerancia fija del 3% como en el paper

while ~isempty(all_poles)
    p = all_poles(1);

    % Tolerancia relativa en real e imaginario (evita división por cero)
    tol_r = tolerance * max(abs(real(p)), 5e-6);
    tol_i = tolerance * max(abs(imag(p)), 5e-6);

    % Encuentra polos cercanos dentro de la tolerancia
    close_idx = abs(real(all_poles) - real(p)) < tol_r & ...
                abs(imag(all_poles) - imag(p)) < tol_i;

    cluster = all_poles(close_idx);

    % Polo representativo: mediana del cluster
    collapsed_poles = [collapsed_poles; median(cluster)];

    % Eliminar el cluster del conjunto
    all_poles(close_idx) = [];
end

collapsed_poles=transpose(collapsed_poles);

    %% Ajuste por subsistema con polos colapsados
    system1 = PCCF1(Yqq, jw, collapsed_poles, weight, pole_tol, red_mode, opts);
    system2 = PCCF1(Yqd, jw, collapsed_poles, weight, pole_tol, red_mode, opts);
    system3 = PCCF1(Yqw, jw, collapsed_poles, weight, pole_tol, red_mode, opts);
    system4 = PCCF1(Ydq, jw, collapsed_poles, weight, pole_tol, red_mode, opts);
    system5 = PCCF1(Ydd, jw, collapsed_poles, weight, pole_tol, red_mode, opts);
    system6 = PCCF1(Ydw, jw, collapsed_poles, weight, pole_tol, red_mode, opts);
    system7 = PCCF1(Ywq, jw, collapsed_poles, weight, pole_tol, red_mode, opts);
    system8 = PCCF1(Ywd, jw, collapsed_poles, weight, pole_tol, red_mode, opts);
    system9 = PCCF1(Yww, jw, collapsed_poles, weight, pole_tol, red_mode, opts);

    %% Construcción de matrices MIMO
    Afull = blkdiag(system1.A, system2.A, system3.A, ...
                    system4.A, system5.A, system6.A, ...
                    system7.A, system8.A, system9.A);

    Bfull = [system1.B, zeros(size(system1.B)), zeros(size(system1.B));
             zeros(size(system2.B)), system2.B, zeros(size(system2.B));
             zeros(size(system3.B)), zeros(size(system3.B)), system3.B;
             system4.B, zeros(size(system4.B)), zeros(size(system4.B));
             zeros(size(system5.B)), system5.B, zeros(size(system5.B));
             zeros(size(system6.B)), zeros(size(system6.B)), system6.B;
             system7.B, zeros(size(system7.B)), zeros(size(system7.B));
             zeros(size(system8.B)), system8.B, zeros(size(system8.B));
             zeros(size(system9.B)), zeros(size(system9.B)), system9.B];

    Cfull = [system1.C, system2.C, system3.C, zeros(size(system4.C)), zeros(size(system5.C)), zeros(size(system6.C)), zeros(size(system7.C)), zeros(size(system8.C)), zeros(size(system9.C));
             zeros(size(system1.C)), zeros(size(system2.C)), zeros(size(system3.C)), system4.C, system5.C, system6.C, zeros(size(system7.C)), zeros(size(system8.C)), zeros(size(system9.C));
             zeros(size(system1.C)), zeros(size(system2.C)), zeros(size(system3.C)), zeros(size(system4.C)), zeros(size(system5.C)), zeros(size(system6.C)),system7.C, system8.C, system9.C];

    Dfull = [system1.D, system2.D, system3.D;
             system4.D, system5.D, system6.D;
             system7.D, system8.D, system9.D];

    %% Modelo estado-espacio
    fitted_ss = ss(Afull, Bfull, Cfull, Dfull, ...
                   'StateName', {}, 'inputname', {'vq_fit', 'vd_fit','w_ref'}, ...
                   'outputname', {'iq_fit', 'id_fit','w_VSC'});
   
    fitted_ss=minreal(fitted_ss);

    if ~isstable(fitted_ss)
        error('The fitted state-space model is unstable. Simulation stopped.');
    end
end
