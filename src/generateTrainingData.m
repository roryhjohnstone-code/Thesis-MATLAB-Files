function [trainedX, trainedY] = generateTrainingData(samplesPerClass)
    if nargin < 1, samplesPerClass = 300; end

    numSensors = 5;
    fs = 100;
    t = linspace(0,10,fs*10);

    classes = ["Normal","Drift","Delamination","Crack","Impact"];

    trainedX = [];               % will become a TABLE on first append
    trainedY = strings(0,1);     % grow as string vector; categorical at the end

    for c = 1:numel(classes)
        for s = 1:samplesPerClass
            baseStrain = 0.5 * sin(2*pi*0.2*t);

            switch classes(c)
                case "Normal"
                    strainData = repmat(baseStrain, numSensors, 1) + 0.01*randn(numSensors, numel(t));
                case "Drift"
                    drift = linspace(0, 0.1, numel(t));
                    strainData = repmat(baseStrain, numSensors, 1) + drift + 0.02*randn(numSensors, numel(t));
                case "Delamination"
                    bump = zeros(size(t)); bump(400:430) = 0.3;
                    strainData = repmat(baseStrain, numSensors, 1) + bump + 0.02*randn(numSensors, numel(t));
                case "Crack"
                    step = zeros(size(t)); step(760:end) = 0.4;
                    strainData = repmat(baseStrain, numSensors, 1) + step + 0.02*randn(numSensors, numel(t));
                case "Impact"
                    spike = zeros(size(t)); spike(447:450) = 0.6;
                    strainData = repmat(baseStrain, numSensors, 1) + spike + 0.02*randn(numSensors, numel(t));
            end

            featTable = extractFBGFeatures(strainData, t);

            % First block seeds the table; subsequent blocks append safely.
            if isempty(trainedX)
                trainedX = featTable;
            else
                trainedX = [trainedX; featTable];  %#ok<AGROW>
            end

            trainedY = [trainedY; repmat(classes(c), numSensors, 1)]; %#ok<AGROW>
        end
    end

    trainedY = categorical(trainedY, classes);
end
