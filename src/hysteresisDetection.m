function detections = hysteresisDetection(signal, threshold, hysteresis_ratio, min_duration)
%   hysteresisDetection - Detects threshold crossings with hysteresis
%   Inputs:
%   signal           - signal vector [1 x N]
%   threshold        - vector [1 x N] or scalar
%   hysteresis_ratio - ratio for lower threshold (e.g., 0.85)
%   min_duration     - minimum sustained length of a detection (samples)
%   Output:
%   detections       - logical vector of detections [1 x N]

    n = length(signal);
    detections = false(1, n);
    isAboveUpper = false;
    consecutive = 0;

    if isscalar(threshold)
        upperThreshold = threshold * ones(1, n);
    else
        upperThreshold = threshold;
    end
    lowerThreshold = upperThreshold * hysteresis_ratio;

    for i = 1:n
        if ~isAboveUpper && signal(i) > upperThreshold(i)
            isAboveUpper = true;
            consecutive = 1;
            fault_start = i;
        elseif isAboveUpper
            if signal(i) > lowerThreshold(i)
                consecutive = consecutive + 1;
            else
                if consecutive >= min_duration
                    detections(fault_start:i-1) = true;
                end
                isAboveUpper = false;
                consecutive = 0;
            end
        end
    end

    if isAboveUpper && consecutive >= min_duration
        detections(n - consecutive + 1:n) = true;
    end
end
