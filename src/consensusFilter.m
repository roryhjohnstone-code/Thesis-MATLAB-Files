function consensusResult = consensusFilter(sensorData)
    % Simple consensus filtering:
    % Returns a binary vector indicating consensus fault detection across sensors per sample

    threshold = 0.5; % Example threshold for "fault" in sensor data

    % Detect samples where any sensor exceeds threshold
    faultIndicators = sensorData > threshold;

    % Consensus = fault present if >= half sensors indicate fault
    numSensors = size(sensorData,2);
    consensusResult = sum(faultIndicators, 2) >= ceil(numSensors/2);
end
