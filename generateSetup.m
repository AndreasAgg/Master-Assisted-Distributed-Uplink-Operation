function [gainOverNoisedB,R,pilotIndex,D,masterAPs,APpositions,UEpositions,distances] = generateSetup(L,K,N,tau_p)
% Generate the simulation setup
% INPUT:
% L               = Number of APs per setup
% K               = Number of UEs in the network
% N               = Number of antennas per AP
% tau_p           = Number of orthogonal pilots
%
% OUTPUT:
% gainOverNoisedB = Matrix with dimension L x K where element (j,k) is 
%                   the channel gain (normalized by the noise variance) 
%                   between AP j and UE k
% R               = Matrix with dimension N x N x L x K
%                   where (:,:,l,k) is the spatial correlation matrix
%                   between AP j and UE k, normalized by noise
% pilotIndex      = Matrix with dimension K x 1 containing the
%                   pilot indices assigned to the UEs
% D               = AP-UE clustering matrix with dimensions L x K where (j,k)
%                   is one if AP j serves UE k and zero otherwise
% APpositions     = Vector of length L with the AP locations, where the real
%                   part is the horizontal position and the imaginary part
%                   is the vertical position
% UEpositions     = Vector of length K with UE positions, measured in the
%                   same way as APpositions
% distances       = Matrix with same dimension as gainOverNoisedB containing
%                   the distances in meter between APs and UEs


%%  Define simulation setup

% Communication bandwidth (Hz)
B = 20e6;

% Noise figure (in dBm)
noiseFigure = 7;

% Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

% Pathloss parameters for the model in (5.42)
alpha = 3.76;

% Standard deviation of shadow fading
sigma_sf = 4;

% Constant term for channel model
constantTerm = -30.5;

% Decorrelation distance of the shadow fading in [12, Eq. (5.43)]
decorr = 9;

% Height difference between an AP and a UE (in meters)
distanceVertical = 10;

% Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; % Half wavelength distance

% Size of the coverage area (as a square with wrap-around)
squareLength = 2000; % meters

% Prepare to save results
gainOverNoisedB = zeros(L,K);
R = zeros(N,N,L,K);
distances = zeros(L,K);
pilotIndex = zeros(K,1);
D = zeros(L,K);
masterAPs = zeros(K,1); % the indices of master AP of each UE k 

% Angular standard deviation in the local scattering model (in radians)
% For more information we refer to [12, Section 2.5]
ASD_varphi = deg2rad(15);  % azimuth angle
ASD_theta = deg2rad(15);   % elevation angle

% Set threshold for when a non-master AP decides to serve a UE
threshold = -40; % dB

%%  Generate setup

% Random AP locations with uniform distribution
APpositions = (rand(L,1) + 1i*rand(L,1)) * squareLength;

% Prepare to compute UE locations
UEpositions = zeros(K,1);


% Compute alternative AP locations by using wrap around
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);

% Prepare to store shadowing correlation matrix
shadowCorrMatrix = sigma_sf^2*ones(K,K);
shadowAPrealizations = zeros(K,L);


% Add UEs
for k = 1:K
    
    % Generate a random UE location in the area
    UEpositions(k) = (rand(1,1) + 1i*rand(1,1)) * squareLength;
    
    % Compute distances assuming that the APs are 10 m above the UEs
    [distanceAPstoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEpositions(k),size(APpositionsWrapped))),[],2);
    distances(:,k) = sqrt(distanceVertical^2+distanceAPstoUE.^2);
    
    % If this is not the first UE
    if k>1
        
        % Compute distances from the new prospective UE to all other UEs
        shortestDistances = zeros(k-1,1);
        
        for i = 1:k-1
            shortestDistances(i) = min(abs(UEpositions(k) - UEpositions(i) + wrapLocations));
        end
        
        % Compute conditional mean and standard deviation necessary to
        % obtain the new shadow fading realizations, when the previous
        % UEs' shadow fading realization have already been generated.
        % This computation is based on Theorem 10.2 in "Fundamentals of
        % Statistical Signal Processing: Estimation Theory" by S. Kay
        newcolumn = sigma_sf^2*2.^(-shortestDistances/decorr);
        term1 = newcolumn'/shadowCorrMatrix(1:k-1,1:k-1);
        meanvalues = term1*shadowAPrealizations(1:k-1,:);
        stdvalue = sqrt(sigma_sf^2 - term1*newcolumn);
        
    else % If this is the first UE
        
        % Add the UE and begin to store shadow fading correlation values
        meanvalues = 0;
        stdvalue = sigma_sf;
        newcolumn = [];
        
    end
    
    % Generate the shadow fading realizations
    shadowing = meanvalues + stdvalue*randn(1,L);
    
    % Compute the channel gain divided by noise power
    gainOverNoisedB(:,k) = constantTerm - 10*alpha*log10(distances(:,k)) + shadowing' - noiseVariancedBm;
    
    
    % Update shadowing correlation matrix and store realizations
    shadowCorrMatrix(1:k-1,k) = newcolumn;
    shadowCorrMatrix(k,1:k-1) = newcolumn';
    shadowAPrealizations(k,:) = shadowing;
    
    % Determine the master AP for UE k by looking for the AP with best
    % channel condition (largest gain)
    [~,master] = max(gainOverNoisedB(:,k));
    D(master,k) = 1;
    masterAPs(k) = master;
    
    % Assign orthogonal pilots to the first tau_p UEs according to
    % [12, Algorithm 4.1]
    if k <= tau_p
        
        pilotIndex(k) = k;
        
    else % Assign pilot for remaining UEs
        
        % Compute received power to the master AP from each pilot
        % according to [12, Algorithm 4.1]
        pilotinterference = zeros(tau_p,1);
        
        for t = 1:tau_p
            
            pilotinterference(t) = sum(db2pow(gainOverNoisedB(master,pilotIndex(1:k-1)==t)));
            
        end
        
        % Find the pilot with the least receiver power according to
        % [12, Algorithm 4.1]
        [~,bestpilot] = min(pilotinterference);
        pilotIndex(k) = bestpilot;
        
    end
    
    
    % Go through all APs
    for j = 1:L
        
        % Compute nominal angle between UE k and AP j
        angletoUE_varphi = angle(UEpositions(k)-APpositionsWrapped(j,whichpos(j))); % azimuth angle
        angletoUE_theta = asin(distanceVertical/distances(j,k));  % elevation angle
        % Generate spatial correlation matrix using the local
        % scattering model in (2.18) and Gaussian angular distribution
        % by scaling the normalized matrices with the channel gain
        R(:,:,j,k) = db2pow(gainOverNoisedB(j,k))*functionRlocalscattering(N,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
    end
    
end


% AP-UE clustering:
% Each AP serves the UE with the strongest channel condition on each of
% the pilots where the AP isn't the master AP, but only if its channel
% is not too weak compared to the master AP

% Go through the APs
for j = 1:L
    
    % Go through the pilot indices
    for t = 1:tau_p

        % Find UEs with the same pilot
        pilotUEs = find(t==pilotIndex);
        
        if sum(D(j,pilotUEs)) == 0 % If the current AP is not a master AP
            
            % Find the UE with pilot t with the best channel
            [gainValue,UEindex] = max(gainOverNoisedB(j,pilotUEs));
            
            % Serve this UE if the channel is at most "threshold" weaker
            % than the master AP's channel
            if gainValue - gainOverNoisedB(masterAPs(pilotUEs(UEindex)),pilotUEs(UEindex)) >= threshold

                D(j,pilotUEs(UEindex)) = 1;

            end
            
        end
        
    end
    
end



