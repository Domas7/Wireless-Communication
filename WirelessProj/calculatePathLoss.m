function pl_dB = calculatePathLoss(distance, frequency, n)
    % Constants
    c = 3e8; % Speed of light in meters per second
    
    % Free-space path loss formula
    pl_dB = 10*n * log10(distance) + 20 * log10(frequency) + 20 * log10(4 * pi / c);
end