%Function checking if stage is moving, from Thorlabs APT sample code 
function r = IsMoving(stage)
    s = stage.GetStatusBits_Bits(0);
    % Read StatusBits returned by GetStatusBits_Bits method and determine if
    % the motor shaft is moving; Return 1 if moving, return 0 if stationary
    r = bitget(abs(s),5)||bitget(abs(s),6);
end
