function [Hhat,H,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p)
% Generate the channel realizations and estimates of these channels for all
% UEs in the entire network. The channels are assumed to follow correlated
% Rayleigh fading and the MMSE estimator is used.
%
% INPUT:
% R                 = Matrix with dimension N x N x L x K where (:,:,j,k) is
%                     the spatial correlation matrix between AP j and UE k,
%                     normalized by noise variance
% nbrOfRealizations = Number of channel realizations (coherence blocks)
% L                 = Number of APs
% K                 = Number of UEs in the network
% N                 = Number of antennas per AP
% tau_p             = Number of pilot signals
% pilotIndex        = Vector containing the pilot assigned to each UE
% p                 = Uplink transmit power per UE (same for everyone)
% 
% OUTPUT:
% Hhat         = Matrix with dimension L*N x nbrOfRealizations x K where
%                (:,n,k) is the estimated collective channel to UE k in
%                channel realization n.
% H            = Matrix with dimension L*N x nbrOfRealizations x K with the
%                true channel realizations. The matrix is organized in the
%                same way as Hhat.
% C            = Matrix with dimension N x N x L x K where (:,:,j,k) is the
%                spatial correlation matrix of the channel estimation error
%                between AP j and UE k, normalized by noise variance

%% Generate channel realizations

% Generate uncorrelated Rayleigh fading channel realizations
H = (randn(L*N,nbrOfRealizations,K)+1i*randn(L*N,nbrOfRealizations,K));


% Go through all channels and apply the spatial correlation matrices
for j = 1:L
    
    for k = 1:K
        
        % Apply correlation to the uncorrelated channel realizations
        Rsqrt = sqrtm(R(:,:,j,k));
        H((j-1)*N+1:j*N,:,k) = sqrt(0.5)*Rsqrt*H((j-1)*N+1:j*N,:,k);
        
    end
    
end


%% Perform channel estimation

% Store identity matrix of size N x N
eyeN = eye(N);

% Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(N,nbrOfRealizations,L,tau_p) + 1i*randn(N,nbrOfRealizations,L,tau_p));

% Prepare to store results
Hhat = zeros(L*N,nbrOfRealizations,K);

% Prepare to store channel error correlation matrices
C = zeros(size(R));

% Go through all APs
for j = 1:L
    
    % Go through all pilots
    for t = 1:tau_p
        
        % Compute processed pilot signal for all UEs that use pilot t
        % according to [12, (4.4)] with an additional scaling factor 
        % \sqrt{tau_p}
        yp = sqrt(p)*tau_p*sum(H((j-1)*N+1:j*N,:,t==pilotIndex),3) + sqrt(tau_p)*Np(:,:,j,t);
        
        % Compute the matrix in [12, (4.6)] that is inverted in the MMSE 
        % channel estimate in (1)
        Psi = (p*tau_p*sum(R(:,:,j,t==pilotIndex),4) + eyeN);
        
        % Go through all UEs that use pilot t
        for k = find(t==pilotIndex)'
            
            % Compute the MMSE channel estimate
            RPsiInv = R(:,:,j,k) / Psi;
            Hhat((j-1)*N+1:j*N,:,k) = sqrt(p)*RPsiInv*yp;


            % Compute the spatial correlation matrix of the estimation
            % error according to [12, (4.9)]
            C(:,:,j,k) = R(:,:,j,k) - p*tau_p*RPsiInv*R(:,:,j,k);

            
        end
        
    end
    
end
