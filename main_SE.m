% This Matlab script generates Figure 1 in the paper:
%
% Andreas Angelou, Pourya Behmandpoor, Marc Moonen "Master-Assisted
% Distributed Uplink Operation for Cell-Free Massive MIMO Networks",
% accepted for publication at ICASSP 2026.

% Empty workspace and close figures
close all;
clear;
clc;
%% Define simulation setup

nbrOfSetups = 20;           % Number of setups to be generated
nbrOfRealizations = 1000;   % Number of coherence blocks
L = 100;                    % Number of APs
N = 4;                      % Number of antennas per AP
K = 40;                     % Number of UEs
p = 100;                    % Power of UEs
tau_c = 200;                % Number of samples per coherence block
tau_p = 10;                 % Number of pilot signals

% ----- For faster simulations use the following -----
% L = 40;
% K = 20;
% nbrOfSetups = 5;
% nbrOfRealizations = 200;
% ----------------------------------------------------
%% Preparations

% Prepare to store spectral efficiency (SE) matricies
SE_C_MMSE_tot = zeros(K, nbrOfSetups); % C-MMSE (centralized)
SE_P_MMSE_tot = zeros(K, nbrOfSetups); % P-MMSE (centralized)
SE_L_MMSE_LSFD_tot = zeros(K, nbrOfSetups); % L-MMSE (distributed)
SE_LP_MMSE_LSFD_tot = zeros(K, nbrOfSetups); % LP-MMSE (distributed)
SE_MADUO_scl_tot = zeros(K, nbrOfSetups); % Proposed scalable (MADUO-scl)
SE_MADUO_tot = zeros(K, nbrOfSetups); % Proposed non-scalable (MADUO)

%% Calculation of SE matrices

% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);

    % Generate one setup with UEs and APs at random locations and assign
    % master APs to the UEs
    [~,R,pilotIndex,D,masterAPs,~,~,~] = generateSetup(L,K,N,tau_p);
    disp(['Generated setup ' num2str(n)]);

    % Generate channel realizations with estimates and estimation error correlation matrices
    [Hhat,H,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    disp(['Channel estimates of setup ' num2str(n)]);
    
    % Compute SE for MADUO, centralized, and distributed operation
    [SE_C_MMSE, SE_P_MMSE, SE_L_MMSE_LSFD, SE_LP_MMSE_LSFD, SE_MADUO_scl, SE_MADUO] = ...
        functionComputeSE_UL(Hhat,H,D,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,masterAPs);

    % Save SE value
    SE_C_MMSE_tot(:,n) = SE_C_MMSE;
    SE_P_MMSE_tot(:,n) = SE_P_MMSE;
    SE_L_MMSE_LSFD_tot(:,n) = SE_L_MMSE_LSFD;
    SE_LP_MMSE_LSFD_tot(:,n) = SE_LP_MMSE_LSFD;
    SE_MADUO_scl_tot(:,n) = SE_MADUO_scl;
    SE_MADUO_tot(:,n) = SE_MADUO;



    %Remove large matrices at the end of analyzing this setup
    clear Hhat H C R;
    
end


%% Plot simulation results

figure
set(gca,'fontsize',12);
hold on; box on; grid on

plot(sort(vec(SE_L_MMSE_LSFD_tot)), linspace(0,1,K*nbrOfSetups), 'LineStyle', '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'DisplayName', 'L-MMSE');
plot(sort(vec(SE_LP_MMSE_LSFD_tot)), linspace(0,1,K*nbrOfSetups), 'LineStyle', '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2.5, 'DisplayName', 'LP-MMSE');
plot(sort(vec(SE_MADUO_tot)), linspace(0,1,K*nbrOfSetups), '-r', 'LineWidth', 1.5, 'DisplayName', 'MADUO');
plot(sort(vec(SE_MADUO_scl_tot)), linspace(0,1,K*nbrOfSetups), '--r', 'LineWidth', 2.5, 'DisplayName', 'MADUO$^{scl}$');
plot(sort(vec(SE_C_MMSE_tot)), linspace(0,1,K*nbrOfSetups), '-k', 'LineWidth', 1.5, 'DisplayName', 'C-MMSE');
plot(sort(vec(SE_P_MMSE_tot)), linspace(0,1,K*nbrOfSetups), '--k', 'LineWidth', 2.5, 'DisplayName', 'P-MMSE');

ttl = ['$$L = ', num2str(L), ', K = ', num2str(K), ...
       ', N = ', num2str(N), '$$'];

title(ttl, 'Interpreter', 'latex')
xlabel('Spectral Efficiency (bits/s/Hz)','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend('Location', 'Best', 'Interpreter','Latex');
