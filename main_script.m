clear all; 
close all; 

addpath('wcmfb_utils');
%% Warped cosine-modulated filter bank (WCMFB) design

% 1. Filter bank parameters:
fs = 16000;         % Sampling frequency
M = 12;              % Number of channels
m = 4;              % Order of polyphase components
N = 2*m*M;          % Order of filter prototype

%1.2 Parameters of design procedure
rho = 0.95;             % Overlap factor of adjacent bands [0...1]
Npt = round(5.5*N);     % Number of points in frequency grid
psi_term = 0.04;        % Constant that determines termination of optimization [0...1].
max_itr = 10;           % The maximum number of optimization steps (13)
gamma   = 1;            % Ripple ratio
%% 2. Selection of warping coefficient
alpha = -(0.1957 - 1.048*((2/pi)*atan(0.07212*(fs/1000))).^(1/2));

%% 3. Computation of subsampling ratios for WCMFB
[Sk] = subsampling_selection(M, -alpha, rho);

%% 4. Compute a rough approximation of the prototype filter
h = fir1(N-1, 1.21/(2*M),'noscale');

%% 5. Filter prototype optimization
name_wfp = ['warped_fp_M' num2str(M) '_h' num2str(N) '.mat'];

h_opt = warped_cmfb_opt_matlab_solver(h, alpha, Sk, 'new', max_itr, Npt, psi_term, gamma);

%% 6. Plotting the main characteristics of the WCMFB
[xr] = warped_cmfb_characteristics(h, h_opt, alpha, M, 'eng', Npt);

%% 7. Phase-compensation filter
C = allpass_char(alpha,N);

figure;
plot((0:length(xr)-1)./length(xr),unwrap(angle(fft(filter(C,1,xr)))),'LineWidth',2,'Color',[0.1 0.6 0.1]);
grid on;hold on;xlim([0 0.5]);
ylabel('$\varphi(\omega), \mathrm{deg}$', 'Interpreter', 'Latex','FontSize',14);
xlabel('$\omega/2\pi$', 'Interpreter', 'Latex','FontSize',14);


figure;
[gd,w]=grpdelay(xr);
ylabel('samples');
xlim([0 0.5]);
plot(w/pi,gd,'LineWidth',2.5,'Color',[0.1 0.6 0.6]); 
grid on; xlim([0 0.5]);
ylabel('Samples','FontSize',12);
xlabel('$\omega / 2\pi$', 'Interpreter', 'Latex','FontSize',14);

figure;
[gd,w]=grpdelay(filter(C,1,xr));
ylabel('samples');
xlim([0 0.5]);
plot(w/pi,gd,'LineWidth',2.5,'Color',[0.1 0.6 0.6]); 
grid on; xlim([0 0.5]);
ylabel('Samples','FontSize',12);
xlabel('$\omega / 2\pi$', 'Interpreter', 'Latex','FontSize',14);

%% 8. Structure with wcmfb parameters
WarpedCMFB.NumberOfChannel = M;
WarpedCMFB.OrderOfFP = N;
WarpedCMFB.OrderOfPPComponents = m;
WarpedCMFB.Alpha = alpha;
WarpedCMFB.Fs = fs;
WarpedCMFB.CoeffOfFP = h_opt;
WarpedCMFB.SubsamplingFactors = Sk;
WarpedCMFB.SubsamplingFactors = ones(1,M);
WarpedCMFB.PhaseCorr = C;
name_param = ['WarpedCMFB_M' num2str(M) '_h' num2str(N) '_fs' num2str(fs)...
              '.mat'];
save(name_param, 'WarpedCMFB');


%% 9. Example of analysis/synthesis

[x, fs] = audioread('test.wav');

name_param = ['WarpedCMFB_M' num2str(M) '_h' num2str(N) '_fs' num2str(fs)...
              '.mat'];
load(name_param);            % Load structure with wcmfb parameters

Npt = length(x);
n=0:Npt-1;

[X] = warped_cmfb_analysis(WarpedCMFB, x);

[Y, sb_pt] = warped_cmfb_decimation(WarpedCMFB, X);

WarpedCMFB.SubbandSamples = sb_pt;
save WarpedCMFB.mat WarpedCMFB;

[X] = warped_cmfb_interpolation(WarpedCMFB, Y);
[x_rec] = warped_cmfb_synthesis(WarpedCMFB, X);

audiowrite('test_reconstructed.wav',x_rec,fs);

figure;
subplot(211);
plot(n,x,n-170,x_rec); xlim([min(n) max(n)]); grid on;
subplot(212);
plot(n,20*log10(abs(fft(x_rec))),'LineWidth',1.2,'color','g'); xlim([min(n) max(n)/2]); grid on;
figure;
subplot(321);
plot(n,X(1,:)); xlim([min(n) max(n)]); grid on;
subplot(322);
plot((n/Npt)*fs,20*log10(abs(fft(X(1,:)))),'LineWidth',1.5); xlim([0 fs/2]); grid on;
subplot(323);
plot(n,X(3,:)); xlim([min(n) max(n)]); grid on;
subplot(324);
plot((n/Npt)*fs,20*log10(abs(fft(X(3,:)))),'LineWidth',1.5); xlim([0 fs/2]); grid on;
subplot(325);
plot(n,X(8,:)); xlim([min(n) max(n)]); grid on;
subplot(326);
plot((n/Npt)*fs,20*log10(abs(fft(X(8,:)))),'LineWidth',1.5); xlim([0 fs/2]); grid on;