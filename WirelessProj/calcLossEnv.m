function power_after = calcLossEnv(in,dist,type,power)
%CALCLOSSENV 
% Calculate the power lost through an environment
%   Detailed explanation goes here
    power_after=power;
    if in==true
        loss=SendSignal(type,dist);
        power_after=power-loss;
    end
end

