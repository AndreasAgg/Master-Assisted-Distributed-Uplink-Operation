% This Matlab script generates Figures 2 and 3 in the paper:
%
% Andreas Angelou, Pourya Behmandpoor, Marc Moonen "Master-Assisted
% Distributed Uplink Operation for Cell-Free Massive MIMO Networks",
% accepted for publication at ICASSP 2026.

%Empty workspace and close figures
close all;
clear;
clc;

%% Define simulation setup

% Number of Monte-Carlo setups
nbrOfSetups = 50;

% Number of APs 
L = 100;

% Number of UEs in the network
K = 20:20:180;

% Number of antennas per AP
N = 4;

% Length of coherence block
tau_c = 200;

% Length of pilot sequences
tau_p = 10;

% Samples for uplink data transmission
tau_u = tau_c - tau_p;

% ---- For faster simulations use the following -----
% nbrOfSetups = 5;
% L = 50;
% K = 20:20:100;
% ---------------------------------------------------

% Prepare to store the AP-UE clustering matrices 
D_tot = zeros(L,max(K),length(K),nbrOfSetups);
masterAPs_tot = zeros(max(K),length(K),nbrOfSetups);
%% Go through all setups
for n = 1:nbrOfSetups

    % Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);

    % Go through all the UE numbers
    for ik=1:length(K)

        % Generate one setup with UEs and APs at random locations
        [~,~,~,D,masterAPs,~,~,~] = generateSetup(L,K(ik),N,tau_p);

        % Save AP-UE clustering matrix
        D_tot(:,1:K(ik),ik,n) = D;
        masterAPs_tot(1:K(ik),ik,n) = masterAPs;
    end

    clear D masterAPs

end

% Compute or prepare to store the total number of complex scalars sent from 
% the APs to the CPU over the fronthaul per coherence block according to
% [12, Table 5.2]
Cent_Front = tau_c*N*L*ones(1,length(K)); % Centralized
Dist_Front_DCC = zeros(1,length(K)); % Distributed

% Prepare to store the total number of complex scalars sent from the
% ASAPs to the MAP over the fronthaul per coherence block according to
% Section 3.3
MADUO_Front = zeros(1,length(K)); % MADUO
MADUO_scl_Front = zeros(1,length(K)); % MADUO-scl

% Compute or prepare to store the total number of complex scalars sent from 
% the APs to the CPU over the fronthaul for each realization channel statistics 
% according to [12, Table 5.2]
opt_LSFD_Front_DCC = zeros(1,length(K)); % opt LSFD 
nopt_LSFD_Front_DCC = zeros(1,length(K)); % n-opt LSFD


% Compute or prepare to store the number of complex multiplications required 
% per coherence block to compute the combining vectors of all UEs for the
% centralized operation according to [12, Table 5.1]
% The common operations will be counted once by modifying the terms in 
% [12, Table 5.1] accordingly.
C_MMSE_Multip = zeros(1,length(K)); % C-MMSE
P_MMSE_Multip = zeros(1,length(K)); % P-MMSE

% Compute or prepare to store the number of complex multiplications required 
% per coherence block to compute the combining vectors of all UEs for the
% dsitributed operation according to [12, Table 5.3]. 
% The common operations will be counted once by modifying the terms in 
% [12, Table 5.3] accordingly.
L_MMSE_Multip = zeros(1,length(K)); % L-MMSE
LP_MMSE_Multip = zeros(1,length(K)); % LP-MMSE


% Prepare to store the number of complex multiplications required 
% per coherence block to compute the combining vectors of all UEs for the
% MADUO 
MADUO_Multip = zeros(1,length(K));
MADUO_scl_Multip = zeros(1,length(K));

% Go through all the UE numbers
for ik = 1:length(K)
    
    % For UE number K(ik), obtain the mean of \sum_{k=1}^K |A_k| 
    % (the mean of \sum_{j=1}^L |U_j|) over all APs and UEs.
    sumUj = mean(D_tot(:,1:K(ik),ik,:),'all')*L*K(ik);
    
    % Compute the fronthaul signaling load per coherence block load for 
    % "Distributed" according to [12, Table 5.2].
    Dist_Front_DCC(ik) = tau_u*sumUj;

    % Compute the fronthaul signaling load per statistics for 
    % "opt LSFD" according to [12, Table 5.2].
    opt_LSFD_Front_DCC(ik) = (3*K(ik)+1)/2*sumUj;

    % Go through all the setups
    for n = 1:nbrOfSetups
        
        % Extract the AP-UE clustering matrix for setup n and UE number K(ik)
        Dn = reshape(D_tot(:,1:K(ik),ik,n), [L, K(ik)]);

        % Number of APs that serve at least one UE
        L_used = length(find(sum(Dn,2)>=1));
        
        % Update the number of complex multiplication for C-MMSE combining
        % according to [12, Table 5.1] (the common operations for the computation
        % of different UEs' combining vectors are counted once)
        C_MMSE_Multip(ik) = C_MMSE_Multip(ik) + ...
            (N*tau_p+N^2)*K(ik)*L_used + ...
            ((N*L_used)^2+N*L_used)/2*K(ik);
        
        % Update the number of complex multiplication for L-MMSE and LP-MMSE
        % combining according to [12, Table 5.3] (the common operations for the
        % computation of different UEs' combining vectors are counted once)
        L_MMSE_Multip(ik) = L_MMSE_Multip(ik) +...
            (N*tau_p+N^2)*K(ik)*L_used + ...
            (N^2+N)/2*K(ik)*L_used+N^2*sumUj+(N^3-N)/3*L_used;
        
        LP_MMSE_Multip(ik) = LP_MMSE_Multip(ik) + ...
            (N*tau_p+N^2)*sumUj + ...
            (N^2+N)/2*sumUj+N^2*sumUj+(N^3-N)/3*L_used;

        % Extract the master APs for setup n and UE number K(ik)
        masterAPs_n = reshape(masterAPs_tot(1:K(ik),ik,n), [K(ik), 1]);


        % Update the number of complex mult. for MADUO and MADUO-scl
        MADUO_Multip(ik) = MADUO_Multip(ik) + ...
            (N*tau_p+N^2)*K(ik)*L_used; % Channel estimation

        MADUO_scl_Multip(ik) = MADUO_scl_Multip(ik) + ...
            (N*tau_p+N^2)*sumUj; % Channel estimation


        % Go through all UEs
        for k = 1:K(ik)
            
            % Find the APs that serve UE k in setup n
            servingAPs = find(Dn(:,k)==1);
            
            % Find the UEs that are served partially by the same set of APs
            % as UE k, i.e., the set in [12, Eq. (5.15)]
            activeUEs = find(sum(Dn(servingAPs,:),1)>=1);

            % Compute the corresping number of UEs and APs for the above sets
            Sk = length(activeUEs);
            Lk = length(servingAPs);

            % Update the fronthaul signaling load per statistics for 
            % "n-opt LSFD" according to [12, Table 5.2].
            nopt_LSFD_Front_DCC(ik) = nopt_LSFD_Front_DCC(ik)+(3*Sk+1)/2*Lk/nbrOfSetups;

            % This is computed to count the number of operations once
            % Number of APs serving at least one active UE
            Bk = length(find(sum(Dn(:,activeUEs),2)>=1));

            % Update the number of complex multiplication for C-MMSE and 
            % P-MMSE, in the centralized operation according to 
            % [12, Table 5.1] (the common operations for the computation
            % of different UEs' combining vectors are counted once)
            C_MMSE_Multip(ik) = C_MMSE_Multip(ik) + ...
                (N*Lk)^2+((N*Lk)^3-N*Lk)/3;
            
            P_MMSE_Multip(ik) = P_MMSE_Multip(ik) + ...
                ((N*tau_p+N^2)*Bk+...
                ((N*Bk)^2+N*Bk)/2+(N*Lk)^2+((N*Lk)^3-N*Lk)/3);

            % Update number of complex mult. for MADUO
            MADUO_Multip(ik) = MADUO_Multip(ik) + ...
                (N+Lk-1)^2 + ...                        % inv(B)*\bar{g} (matrix x vector) -- master AP
                ((N+Lk-1)^3-N-Lk+1)/3 + ...             % inv(B) (matrix inversion) -- master AP
                N*(Lk-1)*(K(ik)-1) + ...                % Off-diagonal blocks of B -- master AP
                ((Lk-1)^2+Lk-1)/2*(K(ik)-1) ...         % G*G' on the down right block of B -- master AP
                + N*K(ik)*(Lk-1) + ...                      % Fuse channels -- ASAPs serving UE k
                (N^2+N)/2*K(ik)*Lk + N^2*(Lk-1) + (N^3-N)/3*(Lk-1) + ...    % RC -- ASAPs serving UE k
                (N^2+N)*(Lk-1);                             % Hermitian product for F_k -- ASAPs serving UE k

            MADUO_scl_Multip(ik) = MADUO_scl_Multip(ik) + ...
                (N+Lk-1)^2 + ...                    % inv(B)*\bar{g} (matrix x vector) -- master AP
                ((N+Lk-1)^3-N-Lk+1)/3 + ...         % inv(B) (matrix inversion) -- master AP
                N^2*(Lk-1) + (N^3-N)/3*(Lk-1) + ... % RC vector (2/3 terms) -- ASAPs serving UE k
                (N^2+N)*(Lk-1);                     % Hermitian product for F_k -- ASAPs serving UE k


            % Master AP (MAP) of UE k
            l = masterAPs_n(k);

            % Going through the serving APs of UE k to update the
            % number of complex mult. for MADUO and MADUO-scl
            for jj=1:length(servingAPs)

                % Number of UEs served by the current AP
                nbr_Uj = nnz(Dn(jj,:));

                % Update the number of computations for MADUO
                MADUO_scl_Multip(ik) = MADUO_scl_Multip(ik) + ...
                    (N^2+N)/2*nbr_Uj; % combining vector (1st term of [12, Table 5.3] for LP-MMSE) -- ASAPs serving UE k


                j = servingAPs(jj); % Current AP

                % Skip if current AP is the MAP
                if j==l
                    continue
                end

                % Find the number of UEs served by both MAP l and ASAP j
                U_jl = find((Dn(j,:) == 1) & (Dn(l,:) == 1));
                nbr_U_jl = length(U_jl);

                % Update the number of computations for MADUO-scl
                MADUO_scl_Multip(ik) = MADUO_scl_Multip(ik) + ...
                    N*(nbr_U_jl-1) + ... % G*H matrix multiplication -- master AP
                    N*nbr_Uj; % fuse channels -- ASAPs serving UE k

                % Go through the rest of serving APs with larger index that AP j
                for qq=jj:length(servingAPs)
                    q = servingAPs(qq); % current AP with larger index than AP j

                    % Find the number of UEs served by both APs j and q
                    U_jq = find((Dn(j,:) == 1) & (Dn(q,:) == 1));
                    nbr_U_jq = length(U_jq);
    
                    % Update the number of computations for MADUO-scl
                    MADUO_scl_Multip(ik) = MADUO_scl_Multip(ik) + ...
                        nbr_U_jq; % G*G -- master AP
                end

            end
        end





        % Fronthaul signaling
        for j=1:L
            % Find the served UEs by AP j
            servedUEs = find(Dn(j,:)==1);
            Uj = length(servedUEs);

            % Number of UEs served by AP j as master AP
            Uj_master = sum(masterAPs_n==j);

            % MADUO and MADUO-scl fronthaul signaling
            MADUO_scl_Front(ik) = MADUO_scl_Front(ik) + (tau_u + Uj + 1)*(Uj - Uj_master)/nbrOfSetups;
            MADUO_Front(ik) = MADUO_Front(ik) + (tau_u + K(ik) + 1)*(Uj - Uj_master)/nbrOfSetups;

        end


    end
end


% Compute the means to plot them
Cent_Front = Cent_Front/L;
Dist_Front_DCC = Dist_Front_DCC/L;
MADUO_Front = MADUO_Front/L;
MADUO_scl_Front = MADUO_scl_Front/L;

C_MMSE_Multip = C_MMSE_Multip/L;
P_MMSE_Multip = P_MMSE_Multip/L;
L_MMSE_Multip = L_MMSE_Multip/L;
LP_MMSE_Multip = LP_MMSE_Multip/L;
MADUO_Multip = MADUO_Multip/L;
MADUO_scl_Multip = MADUO_scl_Multip/L;

%% Plot: Fronthaul load (Figure 2)
figure;
hold on; box on; grid on;
set(gca,'fontsize',12);

semilogy(K, Cent_Front,'kd-','LineWidth',1.5, 'DisplayName', 'Centralized');
semilogy(K, Dist_Front_DCC, 'LineStyle', ':', 'Marker', 's', 'Color', [0.6 0.6 0.6],'LineWidth',1.5, 'DisplayName', 'Distributed');
semilogy(K, MADUO_Front, 'r-x', 'Linewidth', 3, 'DisplayName', 'MADUO')
semilogy(K, MADUO_scl_Front, 'r:x', 'Linewidth', 3, 'DisplayName', 'MADUO$^{scl}$')

ttl = ['Fronthaul: $$L = ', num2str(L), ...
       ', N = ', num2str(N), '$$'];

title(ttl, 'Interpreter', 'latex')
xlabel('Number of UEs ($K$)','Interpreter','Latex');
ylabel('Number of complex scalars (data+pilot)','Interpret','Latex');
legend('Interpreter','Latex','Location','Best');
grid on
ylim([700 inf])
%% Plot: Computational Complexity (Figure 3)
figure;
hold on; box on; grid on;
set(gca,'fontsize',12);

semilogy(K, pow2db(C_MMSE_Multip/nbrOfSetups),'kd-','LineWidth',1.5, 'DisplayName', 'C-MMSE');
semilogy(K, pow2db(P_MMSE_Multip/nbrOfSetups),'kd--','LineWidth',1.5, 'DisplayName', 'P-MMSE');
semilogy(K, pow2db(L_MMSE_Multip/nbrOfSetups),'LineStyle', '-', 'Marker', 's','Color', [0.6 0.6 0.6],'LineWidth',1.5, 'DisplayName', 'L-MMSE');
semilogy(K, pow2db(LP_MMSE_Multip/nbrOfSetups),'LineStyle', '--', 'Marker', 's','Color', [0.6 0.6 0.6],'LineWidth',1.5, 'DisplayName', 'LP-MMSE');
semilogy(K, pow2db(MADUO_Multip/nbrOfSetups),'r-x','LineWidth',3, 'DisplayName', 'MADUO');
semilogy(K, pow2db(MADUO_scl_Multip/nbrOfSetups),'r:x','LineWidth',3, 'DisplayName', 'MADUO$^{scl}$');

ttl = ['Computations: $$L = ', num2str(L), ...
       ', N = ', num2str(N), '$$'];

title(ttl, 'Interpreter', 'latex')
grid on
xlabel('Number of UEs ($K$)','Interpreter','Latex');
ylabel('Number of complex multiplications','Interpreter','Latex');
legend('Interpreter','Latex','Location','Best');

yt = yticks;
yt = yt(mod(yt,10) == 0);
yticks(yt)
yticklabels(arrayfun(@(x) ['10^' num2str(x/10)], yt, 'UniformOutput', false))
