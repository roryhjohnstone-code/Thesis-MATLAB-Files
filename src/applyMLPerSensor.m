function [yPredPerSensor, eventRows] = applyMLPerSensor(M, signals, time)
% M: struct from ensureClassifier() with .predictFcn
% signals: SxN raw strain (NOT pre-normalized)
% time:    1xN

S = size(signals,1); N = size(signals,2);
Fs = round((N-1)/(time(end)-time(1)));

yPredPerSensor = repmat("Normal", S, 1);
eventRows = strings(0,4);

% Grab first contiguous run of *any* event per sensor from detectors in base (if present)
det = struct();
try det = evalin('base','det'); catch, end
hasMasks = isstruct(det) && all(isfield(det,{'impact','crack','delam','drift'}));

for s = 1:S
    if hasMasks
        mk = det.impact(s,:) | det.crack(s,:) | det.delam(s,:) | det.drift(s,:);
    else
        mk = false(1,N);           % no external masks → just use mid-window
    end

    % window center = first detection, else middle of record
    if any(mk), idx = find(mk,1,'first'); else, idx = round(N/2); end

    % build a ~0.6 s window around idx (clamped)
    w = max(4, round(0.6*Fs));
    i1 = max(1, idx - floor(w/2)); i2 = min(N, idx + floor(w/2));

    % features from that window
    feats = ml_extract_features(signals(s, i1:i2), Fs);   % 1xD

    % predict
    yhat = M.predictFcn(feats);
    yPredPerSensor(s) = string(yhat);

    % event row for convenience
    eventRows(end+1,:) = ["Sensor " + s, yPredPerSensor(s), string(time(i1)), string(time(i2))]; %#ok<AGROW>
end
end

