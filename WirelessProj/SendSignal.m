function power_cost = SendSignal(type,distance)
%SENDSIGNAL 
% calculating the power lost through the transmission depending on the type
% of environment
%   Detailed explanation goes here
    if type ==1 %%open env
        shadowing_std_dev = 1; % Standard deviation of log-normal shadowing in dB
        multipath_std_dev = 2; % Standard deviation of multipath fading in dB
        path_loss_exponent=2;
    elseif type ==2 %%forest
        shadowing_std_dev = 4; % Standard deviation of log-normal shadowing in dB
        multipath_std_dev = 2; % Standard deviation of multipath fading in dB
        path_loss_exponent=2.5;
    elseif type ==3  %%mountain
        shadowing_std_dev = 3; % Standard deviation of log-normal shadowing in dB
        multipath_std_dev = 7; % Standard deviation of multipath fading in dB
        path_loss_exponent=3;
    elseif type ==4  %%lake
        shadowing_std_dev = 1; % Standard deviation of log-normal shadowing in dB
        multipath_std_dev = 4; % Standard deviation of multipath fading in dB
        path_loss_exponent=1;
    end
    frequency = 2.4e9; % Signal frequency in Hertz (2.4 GHz for WiFi)
    
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

