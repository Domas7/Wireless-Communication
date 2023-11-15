function power_cost = SendSignal(type,distance)
%SENDSIGNAL Summary of this function goes here
%   Detailed explanation goes here
    frequency = 2.4e9; % Signal frequency in Hertz (2.4 GHz for WiFi)
    shadowing_std_dev = 7; % Standard deviation of log-normal shadowing in dB
    multipath_std_dev = 7; % Standard deviation of multipath fading in dB
    path_loss_exponent=type;
    % Path loss calculation
        
    % Path loss
    pl = calculatePathLoss(distance, frequency, path_loss_exponent);
    
    % Shadowing (log-normal)
    shadowing = normrnd(0, shadowing_std_dev);
    
    % Multipath fading (Rayleigh fading)
    multipath = normrnd(0, multipath_std_dev);
    
    % Received power
    power_cost =  pl + shadowing + multipath;

end

