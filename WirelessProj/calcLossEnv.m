function power_after = calcLossEnv(in,dist,n,power)
%CALCLOSSENV 
% Calculate the power lost through an environment
%   Detailed explanation goes here
    power_after=power;
    if in==true
        loss=SendSignal(n,dist);
        power_after=power-loss;
    end
end

