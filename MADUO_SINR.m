function [SINR_k_MADUO_scl, SINR_k_MADUO] = MADUO_SINR(D,k,l,N,Hhat_n,p,C)
% Compute the SINR given in (9) for the non-scalable and scalable MADUO
% INPUT:
% D                = AP-UE clustering matrix with dimensions L x K where (j,k)
%                    is one if AP j serves UE k and zero otherwise
% k                = index of current UE
%
% l                = AP index of the master AP of UE k
% N                = Number of antennas per AP
% Hhat_n           = Matrix with dimension L*N x K where (:,k) is the 
%                    estimated collective channel to UE k
% p                = Uplink transmit power per UE (same for everyone)
% C                = Matrix with dimension N x N x L x K where (:,:,j,k) is the
%                    spatial correlation matrix of the channel estimation error
%                    between AP j and UE k, normalized by noise variance
%
% OUTPUT:
% SINR_k_MADUO_scl = SINR of UE k for MADUO-scl
% SINR_k_MADUO     = SINR of UE k for MADUO

% Find serving APs of UE k
servingAPs = find(D(:,k)==1);

% Number of UEs in the networks
K = size(D,2);

% Number of APs serving UE k
L_k = length(servingAPs);

% Exclude master AP
ASAPs = setdiff(servingAPs,l);

% Set of all UEs, excluding current UE k
otherUEs = setdiff(1:K,k);

% Prepare to store Ghat_k
Ghat_k_LP_MMSE = zeros(L_k-1, K);
Ghat_k_LP_MMSE_denom = zeros(L_k-1, K-1);
Ghat_k_L_MMSE = zeros(L_k-1, K);
Ghat_k_L_MMSE_denom = zeros(L_k-1, K-1);

% For convinience
eyeN = eye(N);

% Preapre to store F_k
F_k_LP_MMSE = zeros(L_k-1,L_k-1);
F_k_LP_MMSE_denom = zeros(L_k-1,L_k-1);
F_k_L_MMSE = zeros(L_k-1,L_k-1);

% Prepare to store \hat{z}_{l,kk}
zhat_l_kk_LP_MMSE = zeros(N+L_k-1,1);
zhat_l_kk_LP_MMSE(1:N) = Hhat_n((l-1)*N+1:l*N, k);
zhat_l_kk_L_MMSE = zeros(N+L_k-1,1);
zhat_l_kk_L_MMSE(1:N) = Hhat_n((l-1)*N+1:l*N, k);

% Go through all ASAPs
for jj=1:length(ASAPs)
    j = ASAPs(jj); % Current ASAP

    % Extract which UEs are served by ASAP j
    servedUEs = find(D(j,:)==1);

    % Extract channel estimates from all UEs (including k) to ASAP j
    Hhat_j_all = reshape(Hhat_n((j-1)*N+1:j*N,:), [N,K]);


    % Compute LP-MMSE combining according to [12, Eq. (5.39)]
    v_jk_LP_MMSE = p*((p*(Hhat_j_all(:,servedUEs)*Hhat_j_all(:,servedUEs)' + sum(C(:,:,j,servedUEs),4))+eyeN)\Hhat_j_all(:,k));

    % Compute L-MMSE combining according to [12, Eq. (5.29)]
    v_jk_L_MMSE = p*((p*(Hhat_j_all*Hhat_j_all' + sum(C(:,:,j,:),4))+eyeN)\Hhat_j_all(:,k));

    % Assign to Ghat_k defined AFTER (10)
    Hhat_j = reshape(Hhat_n((j-1)*N+1:j*N, otherUEs), [N,K-1]);
    Hhat_j_all_active = Hhat_j_all;
    Hhat_j_all_active(:, D(j,:)==0) = 0;

    Ghat_k_LP_MMSE(jj,:) = v_jk_LP_MMSE'*Hhat_j_all_active;
    Ghat_k_LP_MMSE_denom(jj,:) = v_jk_LP_MMSE'*Hhat_j;
    Ghat_k_L_MMSE(jj,:) = v_jk_L_MMSE'*Hhat_j_all;
    Ghat_k_L_MMSE_denom(jj,:) = v_jk_L_MMSE'*Hhat_j;

    % Assign to F_k defined AFTER (10)
    F_k_LP_MMSE(jj,jj) = v_jk_LP_MMSE'*(p*sum(C(:,:,j,servedUEs),4) + eyeN)*v_jk_LP_MMSE;
    F_k_LP_MMSE_denom(jj,jj) = v_jk_LP_MMSE'*(p*sum(C(:,:,j,:),4) + eyeN)*v_jk_LP_MMSE;
    F_k_L_MMSE(jj,jj) = v_jk_L_MMSE'*(p*sum(C(:,:,j,:),4) + eyeN)*v_jk_L_MMSE;

    % Assign to zhat_l_kk defined AFTER (9)
    zhat_l_kk_LP_MMSE(N+jj) = v_jk_LP_MMSE'*Hhat_n((j-1)*N+1:j*N, k);
    zhat_l_kk_L_MMSE(N+jj) = v_jk_L_MMSE'*Hhat_n((j-1)*N+1:j*N, k);

end

% Collect channel estimates of master AP (excluding UE k)
Hhat_l = reshape(Hhat_n((l-1)*N+1:l*N, otherUEs), [N,K-1]);
Hhat_l_all = reshape(Hhat_n((l-1)*N+1:l*N, :), [N,K]);
Hhat_l_all_active = Hhat_l_all;
Hhat_l_all_active(:, D(l,:)==0) = 0;


% Prepare to store matrix B defined in (10)
B_LP_MMSE = zeros(N+L_k-1, N+L_k-1);
B_L_MMSE = zeros(N+L_k-1, N+L_k-1);
B_LP_MMSE_denom = zeros(N+L_k-1, N+L_k-1);
B_L_MMSE_denom = zeros(N+L_k-1, N+L_k-1);

% Fill matrix B which will be used to calculate the receive combiner in (11)
B_LP_MMSE(1:N,1:N) = p*(Hhat_l_all_active*Hhat_l_all_active' + sum(C(:,:,l,D(l,:)==1),4)) + eyeN;
B_LP_MMSE(1:N, N+1:end) = p*Hhat_l_all_active*Ghat_k_LP_MMSE';
B_LP_MMSE(N+1:end, 1:N) = B_LP_MMSE(1:N, N+1:end)';
B_LP_MMSE(N+1:end, N+1:end) = p*(Ghat_k_LP_MMSE*Ghat_k_LP_MMSE') + F_k_LP_MMSE;

B_L_MMSE(1:N,1:N) = p*(Hhat_l_all*Hhat_l_all' + sum(C(:,:,l,:),4)) + eyeN;
B_L_MMSE(1:N, N+1:end) = p*Hhat_l_all*Ghat_k_L_MMSE';
B_L_MMSE(N+1:end, 1:N) = B_L_MMSE(1:N, N+1:end)';
B_L_MMSE(N+1:end, N+1:end) = p*(Ghat_k_L_MMSE*Ghat_k_L_MMSE') + F_k_L_MMSE;

% Calculating the "true" matrix B which will be used to calculate the SINR in (9)
% The "true" matrix B depends on the channels of ALL UEs in the network
B_LP_MMSE_denom(1:N,1:N) = p*(Hhat_l*Hhat_l' + sum(C(:,:,l,:),4)) + eyeN;
B_LP_MMSE_denom(1:N, N+1:end) = p*Hhat_l*Ghat_k_LP_MMSE_denom';
B_LP_MMSE_denom(N+1:end, 1:N) = B_LP_MMSE_denom(1:N, N+1:end)';
B_LP_MMSE_denom(N+1:end, N+1:end) = p*(Ghat_k_LP_MMSE_denom*Ghat_k_LP_MMSE_denom') + F_k_LP_MMSE_denom;

B_L_MMSE_denom(1:N,1:N) = p*(Hhat_l*Hhat_l' + sum(C(:,:,l,:),4)) + eyeN;
B_L_MMSE_denom(1:N, N+1:end) = p*Hhat_l*Ghat_k_L_MMSE_denom';
B_L_MMSE_denom(N+1:end, 1:N) = B_L_MMSE_denom(1:N, N+1:end)';
B_L_MMSE_denom(N+1:end, N+1:end) = p*(Ghat_k_L_MMSE_denom*Ghat_k_L_MMSE_denom') + F_k_L_MMSE;

% Receive combining at master AP according to (11)
v_k_LP_MMSE = B_LP_MMSE\zhat_l_kk_LP_MMSE;
v_k_L_MMSE = B_L_MMSE\zhat_l_kk_L_MMSE;

% Compute SINR for UE k according to (9)
numerator = p*abs(v_k_LP_MMSE'*zhat_l_kk_LP_MMSE)^2;
denominator = v_k_LP_MMSE'*B_LP_MMSE_denom*v_k_LP_MMSE;
SINR_k_MADUO_scl = numerator / denominator;

numerator = p*abs(v_k_L_MMSE'*zhat_l_kk_L_MMSE)^2;
denominator = v_k_L_MMSE'*B_L_MMSE_denom*v_k_L_MMSE;
SINR_k_MADUO = numerator / denominator;
