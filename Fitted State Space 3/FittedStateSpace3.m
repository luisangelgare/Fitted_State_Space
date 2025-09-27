
% Technical University of Catalonia (UPC)
% Higher Technical School of Industrial Engineering of Barcelona (ETSEIB)
% Centre of Technological Innovation in Static Converters and Drives (CITCEA)
% Doctoral Program in Electrical Engineering
% Developed by: Luis Angel Garcia Reyes, MSc

% State-space methodology for dq control based on PE black-box model

% clear all;
% clc;

%% State-Space identification

Tobs=13; % Observation time
delta_t=25E-6; % Fixed step time
samples=Tobs/delta_t; % Number of samples
time=(0:samples-1)*delta_t; % Time vector 
Delta_f0=1/Tobs; % Frequency step (rad)
f1=(0:samples-1)*Delta_f0; % Frequency vector
jw1=1i*2*pi*f1; % Complex frequency vector
% f0=50; % Base frequency
% w=2*pi*f0; % Base angular frequency (rad/s)
% t_w1=1.8; % Starting time window
% t_w2=2.8; % Ending time window
% Tobs0=t_w2-t_w1; % Total time window
% fs=1/Tobs0; % Sampling frequency for FTT

% load('Yqd_9bus_GFM.mat');

% load('Yqd_GFL_110linear.mat'); % De aqui se extraen las medidiciones de la matriz de admitancia Y y las frecuencias fd0
load('YGFL_3P.mat');

Yqdw_full=[];
% Yqd_full=Yqd0_list{1};
for n=1:length(fd0)
    Yqd_full=[Yqq(n) Yqd(n) Yqw(n)
              Ydq(n) Ydd(n) Ydw(n)
              Ywq(n) Ywd(n) Yww(n)];
    Yqdw_full(:,:,n)=Yqd_full;
end

pole_tol=5e-4; % Numero de iteraciones para los polos
degree=4; % Orden de la seccion 1
fmin=0.1;
fmax=1000;
poles1=-100-1i*linspace(fmin,fmax,degree);
% poles1=-100-1i*logspace(log10(fmin),log10(fmax),degree);
% GFL_VSC = fitss3(Yqdw_full, fd0, poles1, pole_tol, 1);
GFL_VSC = fitss4(Yqdw_full, fd0, [fmin fmax], pole_tol, 1);
% save('TF_VSC_GFL_3P.mat','GFL_VSC','fd0');
% [Ym_Th_gfl0,Ya_Th_gfl0]=bode(SS_GFM0,imag(jw1));
[Ym_Th_gfl0,Ya_Th_gfl0]=bode(GFL_VSC,imag(jw1));
% pole_comparisson(GFL_VSC_ref, GFL_VSC)
% qd0Plot(fd0, jw1, Ym_Th_gfl0, Ya_Th_gfl0, squeeze(Yqd0_full(1,1,:)),...
    % squeeze(Yqd0_full(1,2,:)), squeeze(Yqd0_full(2,1,:)), squeeze(Yqd0_full(2,2,:)));

GFL_VSC_n4sid = MIMOn4sid(Yqdw_full, fd0, delta_t, 3);
isstable(GFL_VSC_n4sid)



clear Yqd_a Ydq_a Ydd_a
for nm=1:length(fd0)
    % if nm<=(length(fd0)-16)
        Ydq_a(nm,:)=angle(Ydq(nm))+2*pi;
    % else
        % Ydq_a(nm,:)=angle(Ydq(nm));
    % end
end

for nm=1:length(fd0)
    % if nm<=(length(fd0)-16)
        Ydd_a(nm,:)=angle(Ydd(nm))+2*pi;
    % else
        % Ydd_a(nm,:)=angle(Ydd(nm));
    % end
end

for nm=1:length(fd0)
    if nm<=24
        Yqd_a(nm,:)=angle(Yqd(nm))+2*pi;
    else
        Yqd_a(nm,:)=angle(Yqd(nm));
    end
end

for nm=1:length(fd0)
    if nm<=63
        Yqw_a(nm,:)=angle(Yqw(nm))+2*pi;
    else
        Yqw_a(nm,:)=angle(Yqw(nm));
    end
end

% Original data
load("Transfer_Multivariable.mat");

[Ym_Th_gfl1,Ya_Th_gfl1]=bode(GFL_VSC_ref,imag(jw1));

% paper_plot_Yqd2(jw1, Ym_Th_gfl0, Ya_Th_gfl0, Ym_Th_gfl1, Ya_Th_gfl1,...
%     fd0, Yqq, Yqd, Ydq, Ydd, Yqd_a, Ydq_a, Ydd_a);
 paper_plot_Yqd3(jw1, Ym_Th_gfl0, Ya_Th_gfl0, Ym_Th_gfl1, Ya_Th_gfl1, fd0,...
                          Yqq, Yqd, Ydq, Ydd, Yqd_a, Ydq_a, Ydd_a, ...
                          Yqw, Ydw, Ywq, Ywd, Yww, Yqw_a);

%  % Supongamos que tienes:
% resp = Ywq;                     % Vector complejo de respuesta
% freq = fd0;        % Frecuencias en Hz
% Ts = delta_t;                         % Sistema continuo
% 
% % Crear objeto idfrd
% G = idfrd(resp, freq, Ts, 'Units', 'Hz');
% 
% % Usar n4sid
% sys = n4sid(G, 2);  % Estima modelo de orden 4
% isstable(sys)
 
 %% Eigenvalues analysis

% save('poles_SVD_case_GFL.mat', 'GFL_VSC_SVD');
% load("poles_SVD_case_GFL.mat");
% [pp,zz]=pzmap(GFL_VSC_ref)
pole_comparisson(GFL_VSC_ref, GFL_VSC)
% pole_comparisson(GFL_VSC_ref, GFL_VSC_n4sid)
% x_limits=[-10000 -9000];
% % x_limits=[-500 0];
% y_limits=[-500 500];
% pole_comp_zoom(GFL_VSC_FULL, GFL_VSC, GFL_VSC_SVD,x_limits, y_limits)

FMODAL(GFL_VSC_ref)

D = eigenproperties_response_luis(GFL_VSC_ref,2,2,[1:1:7])

O = obsv(GFL_VSC_ref.A, GFL_VSC_ref.C);
rank_O = rank(O)
% Asumiendo que ya tienes definido el sistema:
% GFL_VSC_ref = ss(A, B, C, D);

% Supongamos que ya tienes definido tu sistema
sys = GFL_VSC_ref;

% Paso 1: Calcular la matriz de observabilidad
O = obsv(sys.A, sys.C);
rango_O = rank(O);
fprintf('Rango de la matriz de observabilidad: %d\\n', rango_O);

% Paso 2: Obtener transformación de coordenadas
[U, ~, ~] = svd(O);  % SVD para obtener base del subespacio observable
T = U(1:10, 1:rango_O); % Base del subespacio observable

% Paso 3: Transformar el sistema
A_new = T' * sys.A * T;
B_new = T' * sys.B;
C_new = sys.C * T;
D_new = sys.D;

% Paso 4: Crear sistema reducido
sys_reducido = ss(A_new, B_new, C_new, D_new);

% Paso 5: Comparar respuestas
figure;
step(sys, 'r', sys_reducido, 'b--');
legend('Original', 'Reducido');
title('Comparación de respuesta al escalón');

figure;
bode(sys, 'r', sys_reducido, 'b--');
legend('Original', 'Reducido');
title('Comparación de respuesta en frecuencia');


%% Passivity analysis

Pmatrix_fit = squeeze(freqresp(GFL_VSC, fd0, 'Hz')) + pagectranspose(squeeze(freqresp(GFL_VSC, fd0, 'Hz')));
clear lambda_fit
for fn = 1:length(fd0)
    % Calcular la matriz pasiva de Yfull_GFL
    Pmatrix_scan = Yqdw_full(:,:,fn) + pagectranspose(Yqdw_full(:,:,fn));
    lambda_min_scan(fn) = min(eig(Pmatrix_scan));  % Valor propio mínimo para Yfull_GFL
    % Calcular los autovalores del fitted
    lambda_fit(:,fn) = eig(squeeze(Pmatrix_fit(:,:,fn)));
    % Si estamos en la primera frecuencia (fd0(1)), y el autovalor mínimo de la respuesta
    % estimada es positivo, forzar a que sea negativo.
%     if fn == 1 && min(lambda_fit) > 0
%         % Modificar la matriz para que el autovalor mínimo sea positivo
%         Pmatrix_fit(:,:,fn) = Pmatrix_fit(:,:,fn) + eye(size(Pmatrix_fit(:,:,fn))) * abs(min(lambda_fit)) - 1.06; % Hacer que sea positivo
%         lambda_fit(:,fn)=eig(squeeze(Pmatrix_fit(:,:,fn)));
%     end
    lambda_min_fit(fn) = min(lambda_fit(:,fn));  % Valor propio mínimo para la respuesta estimada
end

% Graficación de los resultados
figure
plot(fd0, lambda_min_scan)
hold on
plot(fd0, lambda_min_fit)
plot(fd0, zeros(1, length(fd0)), 'r--')
xlabel('Frequency (Hz)')
ylabel('')
legend('Scanner measurements', 'Fitted reduced order', 'Passive limit')
grid on

%% Time-Domain validation

Tobs=10.0; % Observation time
delta_t=(25E-6); % Fixed step time
samples=Tobs/delta_t; % Number of samples
time=(0:samples-1)*delta_t; % Time vector 
Delta_f0=1/Tobs; % Frequency step (rad)
f1=(0:samples-1)*Delta_f0; % Frequency vector
jw1=1i*2*pi*f1; % Complex frequency vector
f0=50; % Base frequency
w=2*pi*f0; % Base angular frequency (rad/s)
t_w1=1.8; % Starting time window
t_w2=2.8; % Ending time window
Tobs0=t_w2-t_w1; % Total time window
fs=1/Tobs0; % Sampling frequency for FTT
phi=(2/3)*pi; % Phase of 120º between phases
Vq0=0.0; % Steady state q-voltage magnitude
Vd0=0.0; % Steady state d-voltage magnitude
dist_time=100.0; % Disturbance time
Rsource=1E-6; % Equivalent disturbance (Ohm)
droop_act=1; % Droop P-f and Q-v activation time
shortcircuit_time=100*[0.5  0.56]; % Short-circuit interval
freq_step=60; % Time of frequency step (s)
voltage_step=3; % Time of voltage step (s)
load_step=100; % Time of load connection (s)
phase_step=100.56; % Time of phase jump
current_step=100.5;
delta_V=1.01; % Ratio of voltage stepcurrent_step=100.5;
P_step=10; % Time of active power change
Q_step=10; % Time of reactive power change
dist_value_q=0; % Initial disturbance value in q
dist_value_d=0; % Initial disturbance value in d
dist_value_w=0;
dist_value_0=0; % Initial disturbance value in q
dist_value_p=0; % Initial disturbance value in d
dist_value_n=0; % Initial disturbance value in q
Vp0=0; % Steady state q-voltage magnitude
Vn0=0; % Steady state d-voltage magnitude
V00=0; % Steady state d-voltage magnitude
fd=0;
program='Frequency_Scan_VSC_LAGR.slx';

%% General data of the GFL converter

Sbase=2.75E6; % Base power (VSC or System)
Vbase=690; % AC Base voltage (V_L-N)
Vref_DC=3*Vbase; % DC reference voltage (V)
Vpeak=(Vbase/sqrt(3))*sqrt(2); % Peak voltage
Ipeak=Sbase/Vpeak;
fsw=500*f0; % Switching frequency (Hz)
Pref_ss=Sbase*0.5; % Reference active power
Qref_ss=Sbase*0.1; % REference reactive power
delta_w=1.01; % Ratio of frequency change
delta_V=1.01; % Ratio of voltage step
delta_P=0.1; % Ratio of active power change
delta_Q=0; % Ratio of reactive power change
delta_theta=pi/3; % Change in phase jump
delta_I=0.8; % Change in the current injection
Rsc=1E-3; % Fault resistance (Ohm)
Rground=0.01; % Ground fault resistance (Ohm)
R_int_IBGT=1E-3;
R_snubber_IGBT=1E4;
levels=2;

%% Power system parameters

Zbase=(Vbase^2)/Sbase; % Base impedance (ohm) 
Imax=Vbase/Zbase; % Maximun current of VSC
Lbase=Zbase/w; % Base inductance (H)
XR_ratio=3; % X/R ratio
SCR=3; % Short-circuit ratio (SCR)      
Scc=SCR*Sbase; % Short-circuit power (MVA)          
Xcc_th=(Vbase^2)/Scc; % Equivalent reactance (ohm)      
R_th=sqrt(Xcc_th^2/(XR_ratio^2+1)); % Thevenin equivalent resistance (ohm)          
X_th=R_th*XR_ratio; % Thevenin equivalent reactance (ohm)
L_th=X_th/w; % Thevenin equivalent inductance (H)
L_filter=0.25*Lbase; % VSC filter inductance value (H)
R_filter=0.2*Zbase; % VSC filter resistance value (ohm)
% L_filter=Vref_DC/((8*(levels-1)^2)*fsw*0.3*Imax);
% R_filter=(1/3)*sqrt(L_filter);
C_filter=0.15/Zbase/w; % VSC filter capacitance value (F)
R_transf=R_filter/10; % Transformer resistance value (ohm)
L_transf=L_filter/10; % Transformer inductance value (H)
Pload=0.2*Sbase/3; % Load power (W)
Rload=Vbase^2/Pload; % Load resistance (ohm)
Pload_dist=Pload*0.01; % Distribution load power (W)
Rload_dist=Vbase^2/Pload_dist; % Load resistance (ohm)
Lload_dist=Zbase/w; % Load inductance (H)
tau_f=0.001; % Time constant of the ABC-filter
Vindex=Vpeak/Vref_DC; % Voltage index between AC and DC
fc=8*f0;
fcI=6*f0;
% fc=1/(2*pi*R_filter*C_filter);
% fc=1/(2*pi*R_filter*C_filter); % Cut-off frequency of lowpass filter (Hz)
kd_filter=0.707; % Damping factor of the filter
Qfactor=(1/R_filter)*sqrt(L_filter/C_filter); % Quality factor of the RLC filter
Vmag_sw=1; % Switching voltage magnitude (V)
delta_ts=delta_t*10; % 
ps=2*pi*(delta_t/(1/f0));
Vc_ini=-0.9*Vpeak;

%% Transmission line PI model parameters:

line_lenght=1; % Transmission line length (m)
R_line_pu=0.017*line_lenght; % Resistance (pu/m)
X_line_pu=0.092*line_lenght; % Inductive reactance (pu/m)
B_line_pu=0.250*line_lenght; % Capacitive susceptance (pu/m)
R_line=Zbase*R_line_pu; % Line series resistance (ohm)
Xl_line=Zbase*X_line_pu; % Line series reactance (ohm)
L_line=Xl_line/(w*1000); % Line series inductance (H)
Xc_line=(B_line_pu/Zbase); % Line shunt reactance (ohm)
C_line=(Xc_line/w)/2; % Line shunt capacitance (F)

%% Voltage Sourced Converter (VSC) control parameters

% Phase Locked Loop (PLL) (Inner loop)

tau_PLL=0.01; %% PLL time constant (s)
xi_PLL=0.707; % Damping ratio
omega_PLL0=4/(tau_PLL*xi_PLL); % Electrical angular velocity
Kp_PLL=(2*xi_PLL*omega_PLL0)/Vpeak; % Proportional gain of PI
% Kp_PLL=(4*xi_PLL^2)/(Vpeak*tau_PLL) % Proportional gain of the PI control
ts_PLL=(2*xi_PLL)/omega_PLL0; % Time constant 
Ki_PLL=Kp_PLL/ts_PLL; % Intrgral gain of the PI control

% Current Control Loop (Inner loop)

tau_CC=0.1E-3; % Current loop time constant (s)
Kp_CC=L_filter/tau_CC; % Proportional gain of the PI control
Ki_CC=R_filter/tau_CC; % Intrgral gain of the PI control

% PQ Control Loop (Outer loop)

tau_P=0.1; % Active power loop time constant (s)
Kp_P=(1E-3*tau_CC)/tau_P; % Proportial gain of the PI control
Ki_P=1*(1E-3)/tau_P; % Intrgral gain of the PI control
tau_Q=0.1; % Reactive power loop time constant (s)
Kp_Q=1E-3*tau_CC/tau_Q; % Proportial gain of the PI control
Ki_Q=(1*1E-3)/tau_Q; % Intrgral gain of the PI control

% Kp_P=1/tau_P; % Proportial gain of the PI control
% Ki_P=tau_CC/tau_P; % Intrgral gain of the PI control
% tau_Q=tau_P; % Reactive power loop time constant (s)
% Kp_Q=1/tau_Q; % Proportial gain of the PI control
% Ki_Q=tau_CC/tau_Q; % Intrgral gain of the PI control

% P-f and Q-V Droop Loops (Outers loop)

tau_P_droop=0.005; % Cut-off frequency (1/f0) for P filter
Kd_P=20*(Sbase*1E-6)/w; % Proportional P-f gain
tau_Q_droop=0.005; % Cut-off frequency (1/f0) for Q filter
Kd_Q=0.3*50*Sbase*1E-6/Vpeak; % Proportional Q-V gain

out=sim(program);
x0=find(out.tout>=(voltage_step-0.1),1); % Select the linearization point
Vpoc_q0=out.V_POC_q0(x0); % q-grid voltage at the linearization point x0
Vpoc_d0=out.V_POC_d0(x0); % d-grid voltage at the linearization point x0
Vnonlin_poc_q0=out.V_POC_q0; % q-grid voltage at the linearization point x0
Vnonlin_poc_d0=out.V_POC_d0; % d-grid voltage at the linearization point x0
Inon_q=out.I_q_c(x0); % q-VSC current at the linearization point x0
Inon_d=out.I_d_c(x0); % d-VSC current at the linearization point x0
Inonlinear_q=out.I_q_c; % q-VSC current at the linearization point x0
Inonlinear_d=out.I_d_c; % d-VSC current at the linearization point x0
times=out.tout;

Tobs=10.0; % Observation time
delta_t=(25E-6)/4; % Fixed step time
samples=Tobs/delta_t; % Number of samples
time=(0:samples-1)*delta_t; % Time vector 
Delta_f0=1/Tobs; % Frequency step (rad)
f1=(0:samples-1)*Delta_f0; % Frequency vector
jw1=1i*2*pi*f1; % Complex frequency vector
f0=50; % Base frequency
w=2*pi*f0; % Base angular frequency (rad/s)
t_w1=1.8; % Starting time window
t_w2=2.8; % Ending time window
Tobs0=t_w2-t_w1; % Total time window
fs=1/Tobs0; % Sampling frequency for FTT
phi=(2/3)*pi; % Phase of 120º between phases
Vq0=0.0; % Steady state q-voltage magnitude
Vd0=0.0; % Steady state d-voltage magnitude
dist_time=100.0; % Disturbance time
Rsource=1E-6; % Equivalent disturbance (Ohm)
droop_act=1; % Droop P-f and Q-v activation time
shortcircuit_time=100*[0.5  0.56]; % Short-circuit interval
freq_step=60; % Time of frequency step (s)
voltage_step=3; % Time of voltage step (s)

out=sim('Linear2_GFL_LAGR.slx');
Ilinear_q=out.Inonlinear(:,1); % q-VSC current at the linearization point x0
Ilinear_d=out.Inonlinear(:,2); % d-VSC current at the linearization point x0

out=sim('Fit_Linear_GFL_LAGR.slx');
Ilinear_fit_q=out.Inonlinear_fit(:,1); % q-VSC current at the linearization point x0
Ilinear_fit_d=out.Inonlinear_fit(:,2); % d-VSC current at the linearization point x0

% Configuración de los límites
tscs = find((out.tout) >= 2, 1); % Índice del inicio de la ventana temporal
low_axis = out.tout(tscs);      % Límite inferior del eje X
up_axis = out.tout(end);        % Límite superior del eje X

% Subplot 2: Magnitud de I^{POC}_{RMS}(s)
figure
plot(times, Inonlinear_q, 'k', 'LineWidth', 3); % Línea principal
hold on;
plot(out.tout, Ilinear_q, 'b--', 'LineWidth', 3); % Respuesta medida (azul)
plot(out.tout, Ilinear_fit_q, 'r-.', 'LineWidth', 2); % Respuesta medida (rojo)
title('$i^{POC}_{q}$(t)', 'FontSize', 18, 'Interpreter', 'latex'); % Usa LaTeX
ylabel('Magnitude (A)');
xlabel('Time (s)');
legend({'EMT simulation', 'Linearized state-space', 'Fitted state-space'}, 'Location', 'northeast', 'Orientation', 'vertical');
xlim([low_axis up_axis]);
% ylim([321 325]);
grid on; grid minor;
ax = gca;
ax.FontSize = 16;
ax.LineWidth = 2;
ax.LabelFontSizeMultiplier = 1.1;

% figure
% plot(times, Inonlinear_d, 'k', 'LineWidth', 3); % Línea principal
% hold on;
% plot(out.tout, Ilinear_d, 'b--', 'LineWidth', 3); % Respuesta medida (azul)
% plot(out.tout, Ilinear_fit_d, 'r-.', 'LineWidth', 2); % Respuesta medida (rojo)
% title('$i^{POC}_{d}$(t)', 'FontSize', 18, 'Interpreter', 'latex'); % Usa LaTeX
% ylabel('Magnitude (A)');
% xlabel('Time (s)');
% legend({'EMT simulation', 'Linearized state-space', 'Fitted state-space'}, 'Location', 'northeast', 'Orientation', 'vertical');
% xlim([low_axis up_axis]);
% % ylim([321 325]);
% grid on; grid minor;
% ax = gca;
% ax.FontSize = 16;
% ax.LineWidth = 2;
% ax.LabelFontSizeMultiplier = 1.1;

%% Interaction assessment methodology

load("Yqd_grid_stable_110f.mat")
load('Yqd_grid_nostable_110f.mat')

clear M_1 M_2 N_1 N_2
delta_t=25E-6;
for n=1:length(fd0)
    Zgrid_stable=inv([Yqq_g(n) Yqd_g(n) 
                      Ydq_g(n) Ydd_g(n)]);
    Zgrid_unstable=inv([Yqq_g_ns(n) Yqd_g_ns(n) 
                        Ydq_g_ns(n) Ydd_g_ns(n)]);
    Y_GFL=[Yqq(n) Yqd(n) 
           Ydq(n) Ydd(n)];
    % Matlab suggested GNC evaluation
    L_1(:,:,n)=Zgrid_stable*Y_GFL;
    L_2(:,:,n)=Zgrid_unstable*Y_GFL;
    % GNC method (if)
    M_1(n)=det(eye(2)+Zgrid_stable*Y_GFL);
    M_2(n)=det(eye(2)+Zgrid_unstable*Y_GFL);
    % Eigenlocus method (only if)
    N_1(:,n)=eig(eye(2)+Zgrid_stable*Y_GFL);
    N_2(:,n)=eig(eye(2)+Zgrid_unstable*Y_GFL);
end

L_1s=idfrd(L_1,fd0,0);
L_1uns=idfrd(L_2,fd0,delta_t);

% Matlab method
figure
nyquist(L_1s)
hold on;
nyquist(L_1uns)
legend('Stable case', 'Unstable case')

% GNC
figure;
plot(real(M_1), imag(M_1), 'r', 'LineWidth', 1.5, 'MarkerSize', 8); 
hold on
plot(real(M_2), imag(M_2), 'b', 'LineWidth', 1.5, 'MarkerSize', 8); 
plot([min(real(M_1))-1, max(real(M_1))+1], [0, 0], 'k--'); % Eje X
plot([0, 0], [min(imag(M_1))-1, max(imag(M_1))+1], 'k--'); % Eje Y
% quiver(real(M_1(1:end-1)), imag(M_1(1:end-1)), diff(real(M_1)), diff(imag(M_1)), 0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
% quiver(real(M_2(1:end-1)), imag(M_2(1:end-1)), diff(real(M_2)), diff(imag(M_2)), 0, 'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
plot(0, 0, 'cx', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Parte Real');
ylabel('Parte Imaginaria');
title('GNC of the two cases');
legend('Stable case', 'Unstable case')
grid on;

% Eigenlocus -STABLE-
figure;
plot(real(N_1(1,:)), imag(N_1(1,:)), 'g', 'LineWidth', 1.5, 'MarkerSize', 8); 
hold on
plot(real(N_1(2,:)), imag(N_1(2,:)), 'm', 'LineWidth', 1.5, 'MarkerSize', 8); 
plot([min(real(N_1(1,:)))-1, max(real(N_1(1,:)))+1], [0, 0], 'k--'); % Eje X
plot([0, 0], [min(imag(N_1(1,:)))-1, max(imag(N_1(1,:)))+1], 'k--'); % Eje Y
% quiver(real(N_1(1,1:end-1)), imag(N_1(1,1:end-1)), diff(real(N_1(1,:))), diff(imag(N_1(1,:))), 0, 'g', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
% quiver(real(N_1(2,1:end-1)), imag(N_1(2,1:end-1)), diff(real(N_1(2,:))), diff(imag(N_1(2,:))), 0, 'm', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
plot(-1, 0, 'cx', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Parte Real');
ylabel('Parte Imaginaria');
title('Eigenlocus stable case');
legend('Lambda 1', 'Lambda 2')
grid on;

% Eigenlocus -UNSTABLE-
figure;
plot(real(N_2(1,:)), imag(N_2(1,:)), 'g', 'LineWidth', 1.5, 'MarkerSize', 8); 
hold on
plot(real(N_2(2,:)), imag(N_2(2,:)), 'm', 'LineWidth', 1.5, 'MarkerSize', 8); 
plot([min(real(N_2(1,:)))-1, max(real(N_2(1,:)))+1], [0, 0], 'k--'); % Eje X
plot([0, 0], [min(imag(N_2(1,:)))-1, max(imag(N_2(1,:)))+1], 'k--'); % Eje Y
% quiver(real(N_2(1,1:end-1)), imag(N_2(1,1:end-1)), diff(real(N_2(1,:))), diff(imag(N_2(1,:))), 0, 'g', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
% quiver(real(N_2(2,1:end-1)), imag(N_2(2,1:end-1)), diff(real(N_2(2,:))), diff(imag(N_2(2,:))), 0, 'm', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
plot(-1, 0, 'cx', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Parte Real');
ylabel('Parte Imaginaria');
title('Eigenlocus unstable case');
legend('Lambda 1', 'Lambda 2')
grid on;