function R = functionRlocalscattering(N,varphi,theta,ASD_varphi,ASD_theta,antennaSpacing)
% Generate the spatial correlation matrix for the local scattering model,
% defined in [12, (2.18)] for Gaussian angular distribution and ULA
% 
% INPUT:
% N              = Number of antennas
% varphi         = Nominal azimuth angle
% theta          = Nominal elevation angle
% ASD_varphi     = Angular standard deviation around the nominal azimuth angle
%                  (measured in radians)
% ASD_theta      = Angular standard deviation around the nominal elevation angle
%                  (measured in radians)
% antennaSpacing = (Optional) Spacing between antennas (in wavelengths)
% 
% OUTPUT:
% R              = N x N spatial correlation matrix


% Set the antenna spacing if not specified by input
if  nargin < 6
    
    % Half a wavelength distance
    antennaSpacing = 1/2;
    
end


% Prepare to compute the first row of the matrix
firstRow  = zeros(N,1);
firstRow(1) = 1;

%% Go through all the columns of the first row
for column = 2:N
    
    % Distance from the first antenna
    distance = antennaSpacing*(column-1);
    
    % Consider different cases depending on the angle spread
    if (ASD_theta>0) && (ASD_varphi>0)
        
        % Define integrand of [12, (2.18)]
        F = @(delta,epsilon)exp(1i*2*pi*distance*sin(varphi+delta).*cos(theta+epsilon)).*...
            exp(-delta.^2/(2*ASD_varphi^2))/(sqrt(2*pi)*ASD_varphi).*...
            exp(-epsilon.^2/(2*ASD_theta^2))/(sqrt(2*pi)*ASD_theta);
        
        % Compute the integral in [12, (2.18)] by including 20 standard deviations
        firstRow(column) = integral2(F,-20*ASD_varphi,20*ASD_varphi,-20*ASD_theta,20*ASD_theta);
        
    elseif ASD_varphi>0
        
        % Define integrand of [12, (2.18)]
        F = @(delta)exp(1i*2*pi*distance*sin(varphi+delta).*cos(theta)).*...
            exp(-delta.^2/(2*ASD_varphi^2))/(sqrt(2*pi)*ASD_varphi);
        
        % Compute the integral in [12, (2.18)] by including 20 standard deviations
        firstRow(column) = integral(F,-20*ASD_varphi,20*ASD_varphi);
        
    elseif ASD_theta>0
        
        % Define integrand of [12, (2.18)]
        F = @(epsilon)exp(1i*2*pi*distance*sin(varphi).*cos(theta+epsilon)).*...
            exp(-epsilon.^2/(2*ASD_theta^2))/(sqrt(2*pi)*ASD_theta);
        
        % Compute the integral in [12, (2.18)] by including 20 standard deviations
        firstRow(column) = integral(F,-20*ASD_theta,20*ASD_theta);
        
    else
        
        firstRow(column) = exp(1i*2*pi*distance*sin(varphi).*cos(theta));
        
    end
    
end

% Compute the spatial correlation matrix by utilizing the Toeplitz structure
% and scale it properly
R = toeplitz(firstRow);
R = R*(N/trace(R));
