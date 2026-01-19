function [SE_C_MMSE, SE_P_MMSE, SE_L_MMSE, SE_LP_MMSE, ...
    SE_MADUO_scl, SE_MADUO] = functionComputeSE_UL(Hhat,H,D,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,masterAPs)
% Compute uplink spectral efficiency (SE) for centralized operation,
% distributed operation, and MADUO
%
% INPUT:
% Hhat              = Matrix with dimension L*N  x nbrOfRealizations x K
%                     where (:,n,k) is the estimated collective channel to
%                     UE k in channel realization n.
% H                 = Matrix with dimension L*N  x nbrOfRealizations x K
%                     with the true channel realizations. The matrix is
%                     organized in the same way as Hhat.
% D                 = Matrix with dimension L x K where (j,k) is one if AP j
%                     serves UE k and zero otherwise
% C                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                     the spatial correlation matrix of the channel
%                     estimation error between AP l and UE k in setup n,
%                     normalized by noise 
% tau_c             = Length of coherence block
% tau_p             = Length of pilot signals
% nbrOfRealizations = Number of channel realizations (coherence blocks)
% N                 = Number of antennas per AP
% K                 = Number of UEs in the network
% L                 = Number of APs
% p                 = Uplink transmit power per UE (same for everyone)
% masterAPs         = Vector K x 1 containing the master APs of all UEs
% 
% OUTPUT:
% SE_C_MMSE         = Vector K x 1  with the SE achieved with C-MMSE combining (centralized operation)
% SE_P_MMSE         = Vector K x 1  with the SE achieved with P-MMSE combining (centralized operation)
% SE_L_MMSE         = Vector K x 1  with the SE achieved with L-MMSE combining (distributed operation)
% SE_LP_MMSE        = Vector K x 1  with the SE achieved with LP-MMSE combining (distributed operation)
% SE_MADUO_scl      = Vector K x 1  with the SE achieved in MADUO-scl
% SE_MADUO          = Vector K x 1  with the SE achieved in MADUO

%% Preparations

% Prepare to store SE for one setup
SE_C_MMSE = zeros(K,1);
SE_P_MMSE = zeros(K,1);
SE_L_MMSE = zeros(K,1);
SE_LP_MMSE = zeros(K,1);
SE_MADUO_scl = zeros(K,1);
SE_MADUO = zeros(K,1);

% Prepare to store statistical quantities for the distributed operation
% (L-MMSE and LP-MMSE)
gki_L_MMSE = zeros(K,K,L);
gki2_L_MMSE = zeros(K,K,L);
Fk_L_MMSE = zeros(L,K);

gki_LP_MMSE = zeros(K,K,L);
gki2_LP_MMSE = zeros(K,K,L);
Fk_LP_MMSE = zeros(L,K);

% For convinience
eyeN = eye(N);

% Compute the pre-log factor
prelogFactor = (1-tau_p/tau_c);

% For printing progress
mod_factor = 200;

%% Compute statistical quantities for distributed operation

% Go through all channel realizations (coherence blocks)
for n = 1:nbrOfRealizations

    for j = 1:L
        % Extract channel realizations from all UEs to AP j
        Hallj = reshape(H((j-1)*N+1:j*N,n,:),[N K]);
        
        % Extract channel estimates from all UEs to AP j
        Hhatallj = reshape(Hhat((j-1)*N+1:j*N,n,:),[N K]);
        
        % Extract which UEs are served by AP j
        servedUEs = find(D(j,:)==1);

        % Compute MR combining according to [12, Eq. (5.32)]
        V_MR = Hhatallj(:,servedUEs);
        
        % Compute L-MMSE combining according to [12, Eq. (5.29)]
        V_L_MMSE = p*((p*(Hhatallj*Hhatallj'+sum(C(:,:,j,:),4))+eyeN)\V_MR);
        
        % Compute LP-MMSE combining according to [12, Eq. (5.39)]
        V_LP_MMSE = p*((p*(V_MR*V_MR'+sum(C(:,:,j,servedUEs),4))+eyeN)\V_MR);
        
        % Compute the conjugates of the vectors g_{ki} in [12, Eq. (5.23)] 
        % for three combining schemes above for the considered channel 
        % realization
        TemporMatr_L_MMSE = Hallj'*V_L_MMSE;
        TemporMatr_LP_MMSE = Hallj'*V_LP_MMSE;

        % Compute the deterministic values of g_{ki}
        
        % Update the sample mean estimates of the expectations in 
        % [12, Eq. (5.27)]
        Fk_L_MMSE(j,servedUEs) = Fk_L_MMSE(j,servedUEs) + vecnorm(V_L_MMSE).^2/nbrOfRealizations;
        Fk_LP_MMSE(j,servedUEs) = Fk_LP_MMSE(j,servedUEs) + vecnorm(V_LP_MMSE).^2/nbrOfRealizations;
        
        % Update the sample mean estimates of the expectations related to g_{ki} in
        % [12, Eq. (5.23)] to be used in the SE expression of 
        % [12, Theorem 5.4]
        gki_L_MMSE(:,servedUEs,j) = gki_L_MMSE(:,servedUEs,j) + TemporMatr_L_MMSE/nbrOfRealizations;
        gki_LP_MMSE(:,servedUEs,j) = gki_LP_MMSE(:,servedUEs,j) + TemporMatr_LP_MMSE/nbrOfRealizations;
        
        gki2_L_MMSE(:,servedUEs,j) = gki2_L_MMSE(:,servedUEs,j) + abs(TemporMatr_L_MMSE).^2/nbrOfRealizations;
        gki2_LP_MMSE(:,servedUEs,j) = gki2_LP_MMSE(:,servedUEs,j) + abs(TemporMatr_LP_MMSE).^2/nbrOfRealizations;

    end
    
end

% Permute the arrays that consist of the expectations that appear in 
% [12, Theorem 5.4] to compute LSFD vectors and the corresponding SEs
gki_L_MMSE = permute(gki_L_MMSE,[1 3 2]);
gki_LP_MMSE = permute(gki_LP_MMSE,[1 3 2]);
gki2_L_MMSE = permute(gki2_L_MMSE,[1 3 2]);
gki2_LP_MMSE = permute(gki2_LP_MMSE,[1 3 2]);

%% Centralized operation

fprintf('Computing spectral efficiency for centralized operation\n')
for n=1:nbrOfRealizations

    for k=1:K

        % Determine the set of serving APs for UE k
        servingAPs = find(D(:,k)==1); %cell-free setup

        % Compute the number of APs that serve UE k
        L_k = length(servingAPs);

        % Determine which UEs that are served by partially the same set
        % of APs as UE k, i.e., the set in [12, Eq. (5.15)]
        servedUEs = sum(D(servingAPs,:),1)>=1;

        % Extract channel realizations and estimation error correlation
        % matrices for the APs involved in the service of UE k
        Hallj_active = zeros(N*L_k,K);
        Hhatallj_active = zeros(N*L_k,K);
        C_tot_blk = zeros(N*L_k,N*L_k);
        C_tot_blk_partial = zeros(N*L_k,N*L_k);

        for j = 1:L_k
            Hallj_active((j-1)*N+1:j*N,:) = reshape(H((servingAPs(j)-1)*N+1:servingAPs(j)*N,n,:),[N K]);
            Hhatallj_active((j-1)*N+1:j*N,:) = reshape(Hhat((servingAPs(j)-1)*N+1:servingAPs(j)*N,n,:),[N K]);
            C_tot_blk((j-1)*N+1:j*N,(j-1)*N+1:j*N) = sum(C(:,:,servingAPs(j),:),4);
            C_tot_blk_partial((j-1)*N+1:j*N,(j-1)*N+1:j*N) = sum(C(:,:,servingAPs(j),servedUEs),4);
        end


        % ----- Centralized MMSE (C-MMSE) combining ----- 

        % Compute C-MMSE combining according to [12, Eq. (5.11)]
        v_k_MMSE = p*((p*(Hhatallj_active*Hhatallj_active')+p*C_tot_blk+eye(L_k*N))\Hhatallj_active(:,k));

        % Compute numerator and denominator of instantaneous SINR in 
        % [12, Eq. (5.5)] for MMSE combining
        numerator = p*abs(v_k_MMSE'*Hhatallj_active(:,k))^2;
        denominator = p*norm(v_k_MMSE'*Hhatallj_active)^2 + v_k_MMSE'*(p*C_tot_blk+eye(L_k*N))*v_k_MMSE - numerator;

        % Update SE C-MMSE  
        SE_C_MMSE(k) = SE_C_MMSE(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;




        % ----- Partial MMSE (P-MMSE) combining ----- 

        %Compute P-MMSE combining according to [12, Eq. (5.16)]
        v_k_P_MMSE = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p*C_tot_blk_partial+eye(L_k*N))\Hhatallj_active(:,k));
        
        %Compute numerator and denominator of instantaneous SINR in [12, Eq. (5.5)]
        numerator = p*abs(v_k_P_MMSE'*Hhatallj_active(:,k))^2;
        denominator = p*norm(v_k_P_MMSE'*Hhatallj_active)^2 + v_k_P_MMSE'*(p*C_tot_blk+eye(L_k*N))*v_k_P_MMSE - numerator;

        % Update SE P-MMSE
        SE_P_MMSE(k) = SE_P_MMSE(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
       
    end
end

%% Distributed operation

fprintf('Computing spectral efficiency for distributed operation\n')
for k = 1:K
    
    % Determine the set of serving APs for UE k
    servingAPs = find(D(:,k)==1);
    
    % Determine which UEs that are served by partially the same set
    % of APs as UE k, i.e., the set in [12, Eq. (5.15)]
    servedUEs = find(sum(D(servingAPs,:),1)>=1);
    
 
    % Expected value of g_{kk}, scaled by \sqrt{p} for L-MMSE combining
    num_vector = conj(vec(sqrt(p)*gki_L_MMSE(k,servingAPs,k)));
    % Compute the matrix whose inverse is taken in [12, Eq. (5.30)] using 
    % the first- and second-order moments of the entries in the vectors 
    % g_{ki}
    temporMatr = gki_L_MMSE(:,servingAPs,k)'*gki_L_MMSE(:,servingAPs,k);
    denom_matrix = p*(diag(sum(gki2_L_MMSE(:,servingAPs,k),1))...
        +temporMatr-diag(diag(temporMatr)))...
        -num_vector*num_vector'+diag(Fk_L_MMSE(servingAPs,k));
    
    % Compute the opt LSFD according to [12, Eq. (5.30)]
    a_opt_L_MMSE = denom_matrix\num_vector;

    % Compute the SE achieved with opt LSFD and L-MMSE combining according to
    % [12, Eq. (5.25)]
    SE_L_MMSE(k) = prelogFactor*real(log2(1+abs(a_opt_L_MMSE'*num_vector)^2/(a_opt_L_MMSE'*denom_matrix*a_opt_L_MMSE)));



    % Expected value of g_{kk}, scaled by \sqrt{p} for LP-MMSE combining
    num_vector = conj(vec(sqrt(p)*gki_LP_MMSE(k,servingAPs,k)));
    % Compute the denominator matrix to compute SE in [12, Theorem 5.4] 
    % using the first- and second-order moments of the entries in the 
    % vectors g_{ki}
    temporMatr = gki_LP_MMSE(:,servingAPs,k)'*gki_LP_MMSE(:,servingAPs,k);
    denom_matrix = p*(diag(sum(gki2_LP_MMSE(:,servingAPs,k),1))...
        +temporMatr-diag(diag(temporMatr)))...
        -num_vector*num_vector'+diag(Fk_LP_MMSE(servingAPs,k));
    
    % Compute the matrix whose inverse is taken in [12, Eq. (5.41)] using 
    % the first- and second-order moments of the entries in the vectors 
    % g_{ki}
    temporMatr = gki_LP_MMSE(servedUEs,servingAPs,k)'*gki_LP_MMSE(servedUEs,servingAPs,k);
    
    denom_matrix_partial =  p*(diag(sum(gki2_LP_MMSE(servedUEs,servingAPs,k),1))...
        +temporMatr-diag(diag(temporMatr)))...
        -num_vector*num_vector'+diag(Fk_LP_MMSE(servingAPs,k));
    
    
    % Compute the n-opt LSFD according to [12, Eq. (5.41)] for LP-MMSE 
    % combining
    a_nopt_LP_MMSE = denom_matrix_partial\num_vector;
    
    % Compute the SE achieved with n-opt LSFD and LP-MMSE combining 
    % according to [12, Eq. (5.25)]
    SE_LP_MMSE(k) = prelogFactor*real(log2(1+abs(a_nopt_LP_MMSE'*num_vector)^2/(a_nopt_LP_MMSE'*denom_matrix*a_nopt_LP_MMSE)));

end

%% MADUO (proposed operation)

fprintf('Computing spectral efficiency for MADUO\n')
% Go through channel realizations
for n=1:nbrOfRealizations
    
    % Print progress
    if mod(n,mod_factor)==0
        fprintf('Realizations: %d/%d\n',n,nbrOfRealizations)
    end

    % Channel estimates of current channel realizations
    Hhat_n = reshape(Hhat(:,n,:), [L*N K]);

    % Go through all UEs
    for k=1:K

        % Find master AP of UE k
        masterAP = masterAPs(k);

        % Proposed (MADUO and MADUO-scl) SINR
        [SINR_k_MADUO_scl, SINR_k_MADUO]  = ...
            MADUO_SINR(D, k, masterAP, N, Hhat_n, p, C);

        % MADUO SE update
        SE_MADUO_scl(k) = SE_MADUO_scl(k) + prelogFactor*real(log2(1+SINR_k_MADUO_scl))/nbrOfRealizations;
        SE_MADUO(k) = SE_MADUO(k) + prelogFactor*real(log2(1+SINR_k_MADUO))/nbrOfRealizations;


    end

end