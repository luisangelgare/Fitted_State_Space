% Technical University of Catalonia (UPC)
% Higher Technical School of Industrial Engineering of Barcelona (ETSEIB)
% Centre of Technological Innovation in Static Converters and Drives (CITCEA)
% Doctoral Program in Electrical Engineering
% Developed by: Luis Angel Garcia Reyes, MSc

%% State-Space Methodology for dq Control Based on PE Black-Box Model

% This script performs state-space identification for dq control systems derived
% from PE black-box models. It fits a state-space representation for dynamic analysis.

%% Parameters and Initialization

% Define observation parameters
Tobs = 3.0;           % Observation time in seconds
delta_t = 25E-6;      % Fixed step time in seconds
samples = Tobs / delta_t; % Number of samples

% Frequency-related parameters
Delta_f0 = 0.1;  % Frequency step in radians
f1 = (0:samples-1) * Delta_f0; % Frequency vector in Hz
jw1 = 1i * 2 * pi * f1; % Complex frequency vector (s-domain representation)

% Base system frequency
f0 = 50;              % Base frequency in Hz
w = 2 * pi * f0;      % Base angular frequency (rad/s)

% Load pre-computed system data
load('Yqd_GFL_110linear.mat'); % Load the admittance matrix data

%% Assemble Full Admittance Matrix Yqd0_full

Yqd0_full = []; % Initialize the 3D matrix to store frequency responses

% Construct Yqd0_full for all frequency samples
for n = 1:length(fd0)
    % Assemble 2x2 matrix at each frequency step
    Yqd_full = [Yqq(n), Yqd(n);
                Ydq(n), Ydd(n)];
    Yqd0_full(:, :, n) = Yqd_full;
end

%% Vector Fitting and State-Space Model Generation

% Define fitting parameters
pole_tol = 5E-3;      % Tolerance for pole iterations
degree = 6;          % Order of the first section (number of poles)
fmin = 0.1;           % Minimum frequency for poles (rad/s)
fmax = 500;           % Maximum frequency for poles (rad/s)
poles1 = -1 - 1i * linspace(fmin, fmax, degree); % Initial pole guess

% Fit the state-space model using the frequency-domain data
GFL_VSC = fitss2(Yqd0_full, fd0, poles1, pole_tol, 1);

%% Compute Frequency Response and Visualize Results

% Compute the magnitude and phase response using Bode analysis
[Ym_Th_gfl0, Ya_Th_gfl0] = bode(GFL_VSC, imag(jw1));

% Plot and compare the admittance components
qd0Plot(fd0, jw1, Ym_Th_gfl0, Ya_Th_gfl0, ...
    squeeze(Yqd0_full(1,1,:)), ...
    squeeze(Yqd0_full(1,2,:)), ...
    squeeze(Yqd0_full(2,1,:)), ...
    squeeze(Yqd0_full(2,2,:)));

% The qd0Plot function visualizes the fitted and measured admittance components

%% Eigenvalues analysis

[poles_fit,zeros_fit]=pzmap(GFL_VSC);


set(0, 'defaultAxesFontSize', 14);
set(0, 'DefaultLineLineWidth', 1.5);
figure
plot(real(poles_fit),imag(poles_fit),'rx')
xlabel('Real axis')
ylabel('Imaginary axis')
grid on
grid minor
