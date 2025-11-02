function writeRunExports(results, cfg)
% writeRunExports  Save per-run metrics (and log) to /exports
% results: struct produced by runSimulation (fields used below)
% cfg: struct with fields {faultType, injectFaults, useML, useKalman, timestamp, notes}

    % ---------- guard rails ----------
    if nargin < 2, cfg = struct; end
    must = {'faultType','injectFaults','useML','useKalman','timestamp'};
    for k = 1:numel(must)
        if ~isfield(cfg,must{k})
            switch must{k}
                case 'timestamp', cfg.(must{k}) = datetime('now');
                otherwise, cfg.(must{k}) = [];
            end
        end
    end
    if ~isfield(cfg,'notes'), cfg.notes = ""; end

    % ---------- ensure folder ----------
    outDir = fullfile(pwd,'exports');
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end

    % ---------- derive tables ----------
    Ns   = numel(results.yTrue);
    runId = char(string(cfg.timestamp, 'yyyyMMdd_HHmmss'));

    % Per-sensor metrics table
    Sensor      = strcat("Sensor ", string((1:Ns).'));
    TrueLabel   = string(results.yTrue(:));
    PredLabel   = string(results.yPred(:));
    TP          = double(results.TP(:));
    FP          = double(results.FP(:));
    FN          = double(results.FN(:));

    Delay_s = results.DelaySec(:);
    if ~isempty(Delay_s) && isnumeric(Delay_s)
        % prettier: keep NaN as NaN in file
    else
        Delay_s = nan(Ns,1);
    end

    Acc_pct = results.AccuracyPerSensor(:);
    if isempty(Acc_pct)
        Acc_pct = 100*double(TrueLabel == PredLabel);
    end

    % Add run-level metadata columns
    FaultType   = repmat(string(cfg.faultType), Ns, 1);
    InjectFaults= repmat(logical(cfg.injectFaults), Ns, 1);
    UseML       = repmat(logical(cfg.useML), Ns, 1);
    UseKalman   = repmat(logical(cfg.useKalman), Ns, 1);
    RunID       = repmat(string(runId), Ns, 1);
    Timestamp   = repmat(cfg.timestamp, Ns, 1);
    Notes       = repmat(string(cfg.notes), Ns, 1);

    perSensorT = table( ...
        RunID, Timestamp, FaultType, InjectFaults, UseML, UseKalman, Notes, ...
        Sensor, TrueLabel, PredLabel, TP, FP, FN, Delay_s, Acc_pct);

    % Confusion matrix table
    classes = string(results.ClassOrder(:));
    CM      = double(results.ConfusionMatrix);
    if isempty(classes)
        % fallback if not present
        classes = unique([TrueLabel; PredLabel],'stable');
        n = numel(classes); CM = zeros(n); 
        for i = 1:Ns
            r = find(classes==TrueLabel(i),1);
            c = find(classes==PredLabel(i),1);
            if ~isempty(r) && ~isempty(c), CM(r,c)=CM(r,c)+1; end
        end
    end
    CM_T = array2table(CM, 'VariableNames', cellstr(classes));
    CM_T = addvars(CM_T, classes, 'Before', 1, 'NewVariableNames', 'TrueClass');

    % ---------- file names ----------
    base = sprintf('metrics_%s_%s', runId, char(cfg.faultType));
    fileRun   = fullfile(outDir, [base '.csv']);
    fileCM    = fullfile(outDir, [base '_confusion.csv']);
    fileLog   = fullfile(outDir, 'metrics_log.csv');   % rolling log

    % ---------- write files (robust to older MATLAB) ----------
    try
        writetable(perSensorT, fileRun);
    catch ME
        warning('Export:PerRunWrite','Failed to write %s: %s', fileRun, ME.message);
    end

    try
        writetable(CM_T, fileCM);
    catch ME
        warning('Export:CMWrite','Failed to write %s: %s', fileCM, ME.message);
    end

    % Append to rolling log
    try
        if exist(fileLog,'file')
            % R2020b+: writetable(...,'WriteMode','append')
            try
                writetable(perSensorT, fileLog, 'WriteMode','append');
            catch
                % older MATLAB: manual append
                oldT = readtable(fileLog);
                newT = [oldT; perSensorT]; 
                writetable(newT, fileLog);
            end
        else
            writetable(perSensorT, fileLog);
        end
    catch ME
        warning('Export:LogAppend','Failed to append to %s: %s', fileLog, ME.message);
    end

    % Nice console ping
    fprintf('📦 Exported run to:\n  - %s\n  - %s\n  - %s (log)\n', fileRun, fileCM, fileLog);
end
