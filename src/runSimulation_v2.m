
function runSimulation_v2(injectFaults, useML, useKalman, faultType, ...
                          axSensors, axConsensus, ax3D, animateFlag, exportFlag, fig) 

% RUNSIMULATION_V2
% GUI-safe simulation + detectors + (optional) ML + tables + labels + exports.

% -------------------- defaults / guards --------------------
if nargin < 1 || isempty(injectFaults), injectFaults = true; end
if nargin < 2 || isempty(useML),        useML        = true; end
if nargin < 3 || isempty(useKalman),    useKalman    = true; end 
if nargin < 4 || isempty(faultType),    faultType    = "Drift"; end
if nargin < 5 || isempty(axSensors),    axSensors    = gobjects(1,5); end
% axConsensus guard
if nargin < 6 || isempty(axConsensus) || ~all(isgraphics(axConsensus,'axes'))
    if nargin >= 10 && isgraphics(fig,'figure')
        axConsensus = axes('Parent',fig,'Visible','on');
    else
        tmp = figure('Visible','off'); axConsensus = axes('Parent',tmp,'Visible','on');
    end
end

% ax3D guard
if nargin < 7 || isempty(ax3D) || ~all(isgraphics(ax3D,'axes'))
    if nargin >= 10 && isgraphics(fig,'figure')
        ax3D = axes('Parent',fig,'Visible','on');
    else
        tmp = figure('Visible','off'); ax3D = axes('Parent',tmp,'Visible','on');
    end
end

if nargin < 8, animateFlag = false; end 
if nargin < 9, exportFlag  = false; end
if nargin < 10 || ~isgraphics(fig,'figure')
    fig = gcf; if isempty(fig) || ~isgraphics(fig,'figure'), fig = figure('Visible','on'); end
end
faultType = string(faultType);

% -------------------- 1) TIME & SYNTHETIC DATA --------------------
Fs = 20; T = 10; N = T*Fs;
time = linspace(0,T,N);
numSensors = 5;

rawStrainData = zeros(numSensors,N);
for s = 1:numSensors
    rawStrainData(s,:) = sin(2*pi*0.2*time + s) + 0.05*cos(2*pi*0.05*time);
end

if useKalman
for s = 1:size(rawStrainData, 1)
    rawStrainData(s,:) = applyKalmanFilter(rawStrainData(s,:), 0.05, 1e-4);
end
end
[~, N] = size(rawStrainData); %S used later, suppress warning
time = time(:).';
assert(numel(time) == N, 'time length must match size(rawStrainData, 2)');
% --- SNR (noise) from fig appdata (fallback to root) ---
snrDB = Inf;  % default = no added noise
if nargin >= 10 && isgraphics(fig,'figure')
    v = getappdata(fig,'SNR_DB');
    if ~isempty(v), snrDB = v; end
end
if ~isfinite(snrDB)
    % root fallback (legacy)
    v = getappdata(0,'FBG_SNR_DB');
    if ~isempty(v), snrDB = v; end
end
assignin('base', 'lastSignals', rawStrainData);
assignin('base', 'lastTime', time);
% --- Apply additive white Gaussian noise to achieve the requested SNR ---
% SNR here is defined vs. the current rawStrainData (pre-fault), per sensor.
if isfinite(snrDB)
    % compute signal power per sensor, then add noise with the right variance
    for s = 1:size(rawStrainData,1)
        rawStrainData(s,:) = awgn(rawStrainData(s,:),snrDB, 'measured');   
    end
    fprintf('[Noise Injection] Target SNR = %.1f dB (manual)\n', snrDB);
else
    fprintf('[Noise Injection] No noise added (SNR = ∞)\n');
end

if injectFaults
    % --- amplify injected faults for visibility ---
    switch faultType
        case "Drift"
            rawStrainData(1,:) = rawStrainData(1,:) + linspace(0,1.0,N);
        case "Delamination"
            seg = safeRange(80,100,N);
            if ~isempty(seg)
            rawStrainData(2,seg) = rawStrainData(2,seg) + 1.2*sin(0.3*(1:numel(seg)));
            end
        case "Crack"
            seg = safeRange(150,N,N);
            if ~isempty(seg)
                rawStrainData(3,seg) = rawStrainData(3,seg) + 2.0;
            end

        case "Impact"
             seg = safeRange(93,95,N);
            if ~isempty(seg)
                % make a short pulse of length numel(seg)
                bump = linspace(3,2,max(1,numel(seg)));
                rawStrainData(4,seg) = rawStrainData(4,seg) + bump;
            end
    end
end

% -------------------- 2) PLOT 5 SENSOR PANELS --------------------
for s = 1:min(5,numSensors)
    if ~isgraphics(axSensors(s),'axes')
        axSensors(s) = subplot(3,5,s,'Parent',fig);
    end
    ax = axSensors(s);
    cla(ax); hold(ax,'on'); grid(ax,'on');
    plot(ax, time, rawStrainData(s,:), 'b-', 'LineWidth', 1.2, 'DisplayName','Strain');
    title(ax, sprintf('Sensor %d Strain', s));
    xlabel(ax,'Time (s)'); ylabel(ax,'Strain');
   % legend(ax,'Location','best','AutoUpdate','off');
    hold(ax,'off');
end

% -------------------- 3) DETECTORS (signature-flexible) --------------------
driftDetected        = false(numSensors,N);
delaminationDetected = false(numSensors,N);
crackDetected        = false(numSensors,N);
impactDetected       = false(numSensors,N);
spikeDetected        = false(numSensors,N); % optional channel

try
    % Preferred 4-output
    [driftDetected, delaminationDetected, crackDetected, impactDetected] = ...
        computeDetections(rawStrainData, time);
catch
    try
        % 5-output (with spike)
        [driftDetected, spikeDetected, delaminationDetected, crackDetected, impactDetected] = ...
            computeDetections(rawStrainData, time);
    catch MEc
        try
            % struct fallback
            Sd = computeDetections(rawStrainData, time);
            if isstruct(Sd)
                if isfield(Sd,'driftMask'),  driftDetected  = maskFixOuter(Sd.driftMask,  numSensors); end
                if isfield(Sd,'delamMask'),  delaminationDetected = maskFixOuter(Sd.delamMask,numSensors); end
                if isfield(Sd,'crackMask'),  crackDetected  = maskFixOuter(Sd.crackMask,  numSensors); end
                if isfield(Sd,'impactMask'), impactDetected = maskFixOuter(Sd.impactMask, numSensors); end
                if isfield(Sd,'spikeMask'),  spikeDetected  = maskFixOuter(Sd.spikeMask,  numSensors); end
            else
                rethrow(MEc);
            end
        catch MEf
            warning('runSimulation_v2:Detectors', 'Detectors failed: %s', MEf.message);
            driftDetected(:)=false; delaminationDetected(:)=false;
            crackDetected(:)=false; impactDetected(:)=false; spikeDetected(:)=false;
        end
    end
end
% --- enforce mutual exclusivity of masks (priority: Impact > Crack > Delam > Drift)
for s = 1:numSensors
    imp = logical(impactDetected(s,:));
    crk = ~imp & logical(crackDetected(s,:));
    dlm = ~imp & ~crk & logical(delaminationDetected(s,:));
    drf = ~imp & ~crk & ~dlm & logical(driftDetected(s,:));
    impactDetected(s,:)       = imp;
    crackDetected(s,:)        = crk;
    delaminationDetected(s,:) = dlm;
    driftDetected(s,:)        = drf;
end
%fprintf('Det hits (nnz): Drift=%d Delam=%d Crack=%d Impact=%d\n', ...
 %   nnz(driftDetected), nnz(delaminationDetected), nnz(crackDetected), nnz(impactDetected));

% --- DRIFT FALSE-POSITIVE SANITIZER (after exclusivity) ---
% Require sustained positive slope within the detected Drift span
% and kill drift if any higher-priority mask also fires on that sensor.
minHold = max( round(0.40*Fs), 8 );   % at least 0.4 s sustained
minSlope = 0.05;                      % tighten if drift still overfires

for s = 1:numSensors
    dr = logical(driftDetected(s,:));
    if ~any(dr)
        continue
    end

    % 1) If any higher-priority event exists on this sensor, drop drift.
    if any(impactDetected(s,:)) || any(crackDetected(s,:)) || any(delaminationDetected(s,:))
        driftDetected(s,:) = false;
        continue
    end

    % 2) Keep only long runs and check the linear trend inside the run
    d = diff([false dr false]);
    b = find(d==1); e = find(d==-1)-1;

    drClean = false(1, N);
    for k = 1:numel(b)
        seg = b(k):e(k);
        if numel(seg) < minHold
            continue
        end
        % slope of y vs time in this segment
        p = polyfit(time(seg), rawStrainData(s,seg), 1);
        if p(1) >= minSlope
            drClean(seg) = true;
        end
    end
    driftDetected(s,:) = drClean;
end
% --- IMPACT SANITIZER (short, strong, single burst) ---
for s = 1:numSensors
    mk = logical(impactDetected(s,:));
    if ~any(mk), continue, end

    d = diff([false mk false]); 
    b = find(d==1); e = find(d==-1)-1;

    keep = false(1,N); bestPk = -inf;
    maxBurst = 3;  % ≤3 samples expected for an impact
    for k = 1:numel(b)
        seg = b(k):e(k);
        if numel(seg) > maxBurst, continue, end
        pk = max(abs(rawStrainData(s,seg)));
        if pk > bestPk
            bestPk = pk; keep = false(1,N); keep(seg) = true;
        end
    end
    % If nothing passed the short-burst test, drop impact
    impactDetected(s,:) = keep;
end

% --- CRACK SANITIZER (clear step + sustained offset; not overlapping impact) ---
minStep  = 0.9;                       % require this much jump
minHold  = max(round(0.50*Fs), 8);    % sustain ~0.5 s
halfWin  = round(0.5*Fs);             % window for pre/post means

for s = 1:numSensors
    mk = logical(crackDetected(s,:));
    if ~any(mk), continue, end

    d = diff([false mk false]); 
    b = find(d==1); e = find(d==-1)-1;

    keep = false(1,N);
    for k = 1:numel(b)
        idx = b(k);                        % first sample of crack run
        pre  = max(1,idx-halfWin):max(1,idx-1);
        post = idx:min(N,idx+halfWin);
        if numel(pre) < round(0.2*Fs) || numel(post) < round(0.2*Fs)
            continue
        end
        stepAmp = mean(rawStrainData(s,post)) - mean(rawStrainData(s,pre));
        if stepAmp >= minStep
            seg = b(k):min(N, b(k)+minHold-1);
            keep(seg) = true;
        end
    end
    % suppress crack if it overlaps impact (priority Impact > Crack)
    keep = keep & ~impactDetected(s,:);
    crackDetected(s,:) = keep;
end


% after detectors, before the ML block
det = struct('impact',impactDetected,'crack',crackDetected, ...
             'delam',delaminationDetected,'drift',driftDetected);
assignin('base','det',det);
% --- keep only the strongest sensor per class to avoid "everything = Drift"
impactDetected       = keepStrongestSensorPerClass(impactDetected);
crackDetected        = keepStrongestSensorPerClass(crackDetected);
delaminationDetected = keepStrongestSensorPerClass(delaminationDetected);
driftDetected        = keepStrongestSensorPerClass(driftDetected);

% -------------------- 4) MACHINE LEARNING  --------------------
predictedLabels = strings(numSensors,1);
classNames = ["Normal","Drift","Delamination","Crack","Impact"]; %#ok<NASGU>

if useML
    M = ensureClassifier();   % expects struct with fields: mdl, mu, sigma, classes, thresholds
    try
        % If ensureClassifier returned just a model, try to load struct M from MAT
        if ~isstruct(M) || ~isfield(M,'mdl')
            Sload = load('classifier_ML.mat');  % expects variable M
            if isfield(Sload,'M'), M = Sload.M; end
        end

        clf    = M.mdl;                   % ClassificationEnsemble
        muTr   = M.mu;                    % 1xD training mean
        sgTr   = M.sigma;                 % 1xD training std
        fnames = ml_feature_names();      % must match training extractor

        % ---- features per sensor (same extractor as training) ----
        Fs = round((N-1)/(time(end)-time(1)));   % Hz
        X  = zeros(numSensors, numel(fnames));
        for s = 1:numSensors
            X(s,:) = ml_extract_features(rawStrainData(s,:), Fs);   % 1xD
        end

        % ---- z-score with training stats ----
        sgTr(sgTr==0) = 1;
        Z = (X - muTr) ./ sgTr;

        % ---- predict posteriors ----
        [lab, score] = predict(clf, Z);          % score: S x K
        lab = string(lab);
        modelClasses = string(labelsToString(clf.ClassNames));  % 1xK
        predictedLabels = lab(:);

        % ---- thresholds from training (fallback 0.5 if missing) ----
        thrMap = buildThresholdMap(M, modelClasses);

        % ---- decide with per-class threshold + margin vs 2nd-best ----
        margin = 0.08;   % tune 0.05–0.15 if you like
        predictedLabels = repmat("Normal", numSensors, 1);
        for i = 1:numSensors
            [pSorted, idxSorted] = sort(score(i,:), 'descend');
            pmax  = pSorted(1);
            imax  = idxSorted(1);
            p2    = pSorted(min(2,numel(pSorted)));   % second-best (or same if K=1)
            cls   = modelClasses(imax);
            t     = thrMap(char(cls));                % numeric threshold for that class

            if pmax >= t && (pmax - p2) >= margin
                predictedLabels(i) = cls;
            else
                predictedLabels(i) = "Normal";
            end
        end
% ---------- ROC/PR score logging (for batch scripts) ----------
% Build a per-sensor yTrue (same logic you use later)
yTrue_tmp = repmat("Normal", numSensors, 1);
switch string(faultType)
    case "Drift",        yTrue_tmp(1) = "Drift";
    case "Delamination", yTrue_tmp(2) = "Delamination";
    case "Crack",        yTrue_tmp(3) = "Crack";
    case "Impact",       yTrue_tmp(4) = "Impact";
end

% Stash ML outputs onto the figure so the batch runner can harvest them
try
    classNames = string(labelsToString(clf.ClassNames));   % 1xK
    setappdata(fig,'LAST_SCORES', struct( ...
        'scores',      score, ...              % S x K
        'classNames',  classNames(:).', ...    % 1 x K string
        'yTrue',       yTrue_tmp(:) ...        % S x 1 string
    ));
catch
    % non-fatal if logging fails
end
% --------------------------------------------------------------

    catch ME
        warning('runSimulation_v2:ML','ML prediction skipped: %s', ME.message);
        predictedLabels(:) = "Normal";
    end
end

% ---- Gate ML by detector evidence (no promotion) ----
% --- Reconcile per sensor: detector first, then ML, else Normal ---
yPred = strings(numSensors,1);
for s = 1:numSensors
    % detector-derived label by precedence
    detLab = majorityLabelFromMasks( ...
        driftDetected(s,:), delaminationDetected(s,:), ...
        crackDetected(s,:), impactDetected(s,:));

    if detLab ~= "Normal"
        % If any detector fired, trust detector precedence
        yPred(s) = detLab;
    else
        % No detector fired: allow ML if non-Normal, else Normal
        if useML && predictedLabels(s) ~= ""
            yPred(s) = predictedLabels(s);
        else
            yPred(s) = "Normal";
        end
    end
end

% ---- Override ML predictions with detector precedence ----
yPred = overrideWithDetectors(yPred, driftDetected, delaminationDetected, crackDetected, impactDetected);

% -------------------- 5) TRUTH / PREDICTIONS --------------------
yTrue = repmat("Normal", numSensors, 1);
switch faultType
    case "Drift",        yTrue(1) = "Drift";
    case "Delamination", yTrue(2) = "Delamination";
    case "Crack",        yTrue(3) = "Crack";
    case "Impact",       yTrue(4) = "Impact";
end

disp("Post-gate yPred: " + strjoin(cellstr(yPred.'), ", "));
% Log scores for ROC (only when ML ran & we have scores)
try
    maybe_log_scores(score, modelClasses, yTrue, yPred, faultType, fig, useML, useKalman);
catch ME
    warning('runSimulation_v2:ScoresLog','Skipping score log: %s', ME.message);
end

% -------------------- 6) CONSENSUS PLOT (class-coloured markers) --------------------
try
    plotConsensusWithClassMarkers(axConsensus, time, rawStrainData, ...
        driftDetected, delaminationDetected, crackDetected, impactDetected);
catch ME
    warning('runSimulation_v2:Consensus','Consensus plot skipped: %s', ME.message);
end


% -------------------- 3D BLADE VIEW --------------------
try
    if isgraphics(ax3D,'axes')
        dm = maskFixOuter(driftDetected,        numSensors);
        dl = maskFixOuter(delaminationDetected, numSensors);
        ck = maskFixOuter(crackDetected,        numSensors);
        im = maskFixOuter(impactDetected,       numSensors);
        sensorStatus = any(dm | dl | ck | im, 2);

        faultMask = struct();
        sensorSpanPos = linspace(0.15,0.85,numSensors).';

        if any(any(dl,2))
            sIdx = find(any(dl,2), 1, 'first');
            faultMask.delamSpan = [ ...
                max(0, sensorSpanPos(max(1,sIdx-1)) - 0.02), ...
                min(1, sensorSpanPos(min(numSensors,sIdx+1)) + 0.02)];
        end
        if any(any(ck,2))
            sIdx = find(any(ck,2), 1, 'first');
            faultMask.crackSpan = sensorSpanPos(sIdx) + [-0.005 0.005];
        end
        if any(any(im,2))
            sIdx = find(any(im,2), 1, 'first');
            faultMask.impactSpan = sensorSpanPos(sIdx) + [-0.005 0.005];
        end

        steps  = iff(animateFlag, 36, 1);
        frames = getappdata(fig,'BLADE_FRAMES');

        if exist('drawTurbine3D_advanced','file')==2
            frames = drawTurbine3D_advanced( ...
                ax3D, sensorStatus, animateFlag, exportFlag, frames, steps, ...
                rawStrainData, [], sensorSpanPos, faultMask);
            if isgraphics(fig,'figure')
    setappdata(fig,'BLADE_FRAMES', frames);
            end

        else
            cla(ax3D);
            imagesc(ax3D, time, 1:numSensors, rawStrainData);
            axis(ax3D,'xy'); xlabel(ax3D,'Time (s)'); ylabel(ax3D,'Sensor');
            title(ax3D,'Strain Heatmap'); colorbar(ax3D);
        end

        % Only adjust Position if ax3D is NOT inside a tiledlayout
        if isempty(ancestor(ax3D,'matlab.graphics.layout.TiledChartLayout'))
            pos = get(ax3D,'Position');
            pos(2) = max(0.05, pos(2)-0.02);
            pos(4) = min(0.28, pos(4));
            set(ax3D,'Position',pos);
        end
    end
catch E
    warning(E.identifier, '%s', E.message);
end

% -------------------- 8) LABEL MARKERS ON SENSOR PLOTS --------------------
driftM  = maskFixOuter(driftDetected,        numSensors);
delamM  = maskFixOuter(delaminationDetected, numSensors);
crackM  = maskFixOuter(crackDetected,        numSensors);
impactM = maskFixOuter(impactDetected,       numSensors);

for s = 1:numSensors
    if s <= numel(axSensors) && isgraphics(axSensors(s),'axes')
        decorateSensorAxes(axSensors(s), time, rawStrainData(s,:), ...
            driftM(s,:), delamM(s,:), crackM(s,:), impactM(s,:), s);
    end
end

% -------------------- 9) DETECTION TIMES (for delays) --------------------
detectedTimes = cell(numSensors,1);
for s = 1:numSensors
    mk = false(1,N);
    if any(impactM(s,:)),      mk = impactM(s,:);
    elseif any(crackM(s,:)),   mk = crackM(s,:);
    elseif any(delamM(s,:)),   mk = delamM(s,:);
    elseif any(driftM(s,:)),   mk = driftM(s,:);
    end
    if any(mk), detectedTimes{s} = time(find(mk,1,'first')); else, detectedTimes{s} = []; end
end

% -------------------- 10) TABLES --------------------
safePopulateFaultTable( ...
    fig, time, ...
    driftDetected, spikeDetected, delaminationDetected, crackDetected, impactDetected, ...
    predictedLabels);

det = struct('impact',impactDetected,'crack',crackDetected, ...
             'delam',delaminationDetected,'drift',driftDetected);
safeUpdateMetricsTable(fig, [], yPred, time, det, detectedTimes);

% -------------------- 11) simResults (+ assignin) --------------------
order = labelsToString(unique([yTrue; yPred],'stable'))';
n  = numel(order);
CM = zeros(n, n);
for i = 1:numSensors
    r = find(order==yTrue(i),1);
    c = find(order==yPred(i),1);
    if ~isempty(r) && ~isempty(c), CM(r,c) = CM(r,c)+1; end
end
TP = diag(CM); 
FP = sum(CM,1)' - TP; %#ok<NASGU>
FN = sum(CM,2) - TP; %#ok<NASGU>
accPerSensor = double(yTrue==yPred)*100;

simResults = struct();
simResults.yTrue              = yTrue(:);
simResults.yPred              = yPred(:);
simResults.DelaySec           = cellfun(@(v) iff(isempty(v), NaN , v - time(1)), detectedTimes);
simResults.AccuracyPerSensor  = accPerSensor;
simResults.AccuracyOverall    = mean(yTrue==yPred);
simResults.ConfusionMatrix    = CM;
simResults.ClassOrder         = order;
% --- EXTRA FIELDS so downstream exports have everything ---
simResults.FaultType     = faultType;
simResults.Time          = time(:).';
simResults.RawStrain     = rawStrainData;
simResults.DetMasks      = struct('drift',driftDetected,...
                                  'delam',delaminationDetected,...
                                  'crack',crackDetected,...
                                  'impact',impactDetected);
% If ML was executed, stash the raw scores + class order
if exist('score','var') && exist('modelClasses','var')
    simResults.MLScores    = score;                % S x K
    simResults.MLClasses   = string(modelClasses); % 1 x K
else
    simResults.MLScores  = [];
    simResults.MLClasses = [];
end
assignin('base','simResults', simResults);

% -------------------- 12) EXPORTS (optional) ------------------------------
if exportFlag
    try
        if exist('writeRunExports','file')
            cfg = struct('faultType',faultType,'injectFaults',injectFaults,...
                         'useML',useML,'useKalman',useKalman,'timestamp',datetime('now'));
            writeRunExports(simResults, cfg);
        end
    catch ME
        warning('runSimulation_v2:ExportMetrics','Export metrics skipped: %s', ME.message);
    end
    try
        if exist('writeFaultSummaryExport','file')
            % your function signature is (fig, outDir, faultType)
            writeFaultSummaryExport(fig, [], faultType);
        end
    catch ME
        warning('runSimulation_v2:ExportFaults','Fault summary export failed: %s', ME.message);
    end
end
end

% ===================== LOCAL HELPERS =====================
function Mout = maskFixOuter(M, numSensors)
% Force masks to SxN (S=numSensors). Accepts SxN or NxS; else returns best-effort.
    if isempty(M)
        Mout = false(numSensors, 0);
        return
    end
    if size(M,1) == numSensors
        Mout = logical(M);
        return
    end
    if size(M,2) == numSensors
        Mout = logical(M.');
        return
    end
    % Best-effort reshape to SxN
    v = logical(M(:));
    if mod(numel(v), numSensors) == 0
        Mout = reshape(v, numSensors, []);
    else
        % pad with false
        n = ceil(numel(v)/numSensors)*numSensors;
        v(end+1:n) = false;
        Mout = reshape(v, numSensors, []);
    end
end

function v = iff(cond, a, b)
    if cond, v = a; else, v = b; end
end

%function lab = majorityLabelFromMasks(driftM, delamM, crackM, impactM)
% precedence: Impact > Crack > Delamination > Drift > Normal
 %   if any(impactM), lab = "Impact";       return; end
  %  if any(crackM),  lab = "Crack";        return; end
   % if any(delamM),  lab = "Delamination"; return; end
    %if any(driftM),  lab = "Drift";        return; end
    %lab = "Normal";
%end

function T = tocharcell(x)
% Convert to cellstr (char) for uitable; accepts string, char, cellstr, numeric, logical.
    if isempty(x), T = {}; return; end
    if isstring(x),  T = cellstr(x); return; end
    if ischar(x),    T = cellstr(x); return; end
    if iscellstr(x) || isstring(x), T = cellstr(x); return; end
    if isnumeric(x)
        T = cellfun(@(v) sprintf('%g',v), num2cell(x), 'UniformOutput', false); return
    end
    if islogical(x)
        T = cellfun(@(v) sprintf('%d',v), num2cell(x), 'UniformOutput', false); return
    end
    T = {};
end

function h = findTable(fig, whichOne)
% Try to find your existing uitables. If not found, create a minimal one.
% whichOne: 'faults' or 'metrics'
    h = [];
    if ~isgraphics(fig,'figure'), return; end
    Ts = findall(fig,'Type','uitable');
    if ~isempty(Ts)
        pos = arrayfun(@(t)get(t,'Position'), Ts, 'uni',0);
        pos = vertcat(pos{:});
        [~,ord] = sort(pos(:,2),'descend');
        Ts = Ts(ord);
        if strcmpi(whichOne,'faults'),  h = Ts(1); return; end
        if numel(Ts) >= 2 && strcmpi(whichOne,'metrics'), h = Ts(2); return; end
    end
    switch lower(whichOne)
        case 'faults'
            h = uitable('Parent',fig,'Units','normalized','Position',[0.72 0.70 0.25 0.22]);
            set(h,'ColumnName',{'Sensor','TP','FP','FN','Delay (s)','Accuracy (%)'});
        case 'metrics'
            h = uitable('Parent',fig,'Units','normalized','Position',[0.72 0.42 0.25 0.22]);
            set(h,'ColumnName',{'Sensor','Fault Type','Start Time (s)','End Time (s)','Delay (s)'});
    end
end

function safePopulateFaultTable(fig, time, driftM, spikeM, delamM, crackM, impactM, predictedLabels)
    S = size(driftM,1);
    N = size(driftM,2);

    driftM  = logical(fixMaskSize(driftM,  S, N));
    spikeM  = logical(fixMaskSize(spikeM,  S, N));
    delamM  = logical(fixMaskSize(delamM,  S, N));
    crackM  = logical(fixMaskSize(crackM,  S, N));
    impactM = logical(fixMaskSize(impactM, S, N));

    ok = false;
    if exist('populateFaultTable','file') == 2
        try
            populateFaultTable(fig, S, time(:).', driftM, spikeM, delamM, crackM, impactM);
            ok = true;
        catch
            try
                detMasks = struct('driftMask',driftM,'spikeMask',spikeM, ...
                                  'delamMask',delamM,'crackMask',crackM,'impactMask',impactM);
                lab = tocharcell(predictedLabels);
                populateFaultTable(fig, time(:).', detMasks, lab);
                ok = true;
            catch
                try
                    populateFaultTable(fig, S, time(:).', driftM, delamM, crackM, impactM, tocharcell(predictedLabels));
                    ok = true;
                catch
                    ok = false;
                end
            end
        end
    end
    if ok, return; end

    labOrder = {'Impact','Crack','Delamination','Drift'};
    rowSensor = (1:S).';
    TP = zeros(S,1); FP = zeros(S,1); FN = zeros(S,1); Delay = nan(S,1);

    for s = 1:S
        mk = impactM(s,:) | crackM(s,:) | delamM(s,:) | driftM(s,:);
        if any(mk)
            Delay(s) = time(find(mk,1,'first')) - time(1);
        end
        p = 'Normal';
        if any(impactM(s,:)),      p = 'Impact';
        elseif any(crackM(s,:)),   p = 'Crack';
        elseif any(delamM(s,:)),   p = 'Delamination';
        elseif any(driftM(s,:)),   p = 'Drift';
        end
        if ~isempty(predictedLabels)
            try
                mlp = string(predictedLabels{s});
                if mlp ~= "", p = char(mlp); end
            catch
            end
        end
        t = 'Normal';
        for k = 1:numel(labOrder)
            switch labOrder{k}
                case 'Impact',      if any(impactM(s,:)),   t = 'Impact';      break; end
                case 'Crack',       if any(crackM(s,:)),    t = 'Crack';       break; end
                case 'Delamination',if any(delamM(s,:)),    t = 'Delamination';break; end
                case 'Drift',       if any(driftM(s,:)),    t = 'Drift';       break; end
            end
        end
        if strcmpi(p,t) && ~strcmpi(t,'Normal')
            TP(s) = 1;
        elseif ~strcmpi(p,t) && ~strcmpi(p,'Normal')
            FP(s) = 1;
        elseif ~strcmpi(p,t) && strcmpi(p,'Normal') && ~strcmpi(t,'Normal')
            FN(s) = 1;
        end
    end

    Acc = 100*(TP>=1);
    data = [ tocharcell("Sensor "+string(rowSensor)), num2cell(TP), num2cell(FP), ...
             num2cell(FN), num2cell(Delay), num2cell(Acc) ];

    ht = findTable(fig,'faults');
    set(ht,'Data',data, ...
        'ColumnName',{'Sensor','TP','FP','FN','Delay (s)','Accuracy (%)'});
end

function safeUpdateMetricsTable(fig, ~, yPred, time, det, ~)
    S = numel(yPred);
    startT = nan(S,1); endT = nan(S,1); delayT = nan(S,1);
    list = {'impact','crack','delam','drift'};
    for s = 1:S
        mk = [];
        for k = 1:numel(list)
            mk = det.(list{k})(s,:);
            if any(mk), break; end
        end
        if any(mk)
            idxStart = find(mk,1,'first'); idxEnd = find(mk,1,'last');
            startT(s) = time(idxStart); endT(s) = time(idxEnd);
            delayT(s) = startT(s) - time(1);
        end
    end
    sensors = "Sensor "+string((1:S).');
    data = [ tocharcell(sensors), tocharcell(yPred), num2cell(startT), num2cell(endT), num2cell(delayT) ];

    ht = findTable(fig,'metrics');
    set(ht,'Data',data, ...
        'ColumnName',{'Sensor','Fault Type','Start Time (s)','End Time (s)','Delay (s)'});
end

function M = fixMaskSize(M,S,N)
    if isempty(M), M = false(S,N); return; end
    M = logical(M);
    if isequal(size(M),[S,N]), return; end
    if isequal(size(M),[N,S]), M = M.'; return; end
    M = M(:).';
    if numel(M) < S*N
        M = [M false(1,S*N-numel(M))];
    end
    M = reshape(M,[S,N]);
end

function decorateSensorAxes(ax, t, y, driftM, delamM, crackM, impactM, sIdx)
% Wrapper: if your plotSensorStrain exists, call it; else draw simple markers.
    if exist('plotSensorStrain','file')==2
        try
            plotSensorStrain(ax, t, y, driftM, [], delamM, crackM, impactM, sIdx);
            return
        catch
            % fall back to simple overlay below
        end
    end
    axes(ax); hold(ax,'on'); grid(ax,'on');
    delete(findobj(ax, 'Type', 'Legend'));
    plot(ax, t, y, 'b-', 'LineWidth', 1.2);
    if any(driftM),  plot(ax, t(driftM),  y(driftM),  'r.',  'MarkerSize',12, 'DisplayName','Drift'); end
    if any(delamM),  plot(ax, t(delamM),  y(delamM),  'gs',  'MarkerSize',5,  'MarkerFaceColor',[0.3 0.8 0.3], 'LineStyle','none', 'DisplayName','Delamination'); end
    if any(crackM)
        d = diff([false crackM false]); b = find(d==1,1,'first');
        if ~isempty(b), plot(ax, t(b), y(b), 'kd', 'MarkerSize',7, 'MarkerFaceColor','k', 'DisplayName','Crack'); end
    end
    if any(impactM), plot(ax, t(impactM), y(impactM), 'p', 'MarkerSize',9, 'LineStyle','none', 'DisplayName','Impact'); end
    title(ax, sprintf('Sensor %d Strain', sIdx));
    xlabel(ax,'Time (s)'); ylabel(ax,'Strain');
    hold(ax,'off');
end
function M = keepStrongestSensorPerClass(M)
% Keep only the single sensor with the longest run (most 1's). Zero others.
    if isempty(M) || ~any(M(:)), return; end
    len = sum(M, 2);              % run length per sensor
    [~, sBest] = max(len);
    keepRow = false(size(M,1),1);
    keepRow(sBest) = true;
    M(~keepRow, :) = false;
end
function thrMap = buildThresholdMap(M, modelClasses)
% Build containers.Map of per-class thresholds aligned to clf.ClassNames.
% Falls back to 0.5 for any missing class.
    keys = cellstr(modelClasses);
    vals = num2cell(0.5*ones(1,numel(keys)));  % default

    % If training stored thresholds, use them
    try
        if isfield(M,'thresholds') && isfield(M,'classes')
            trKeys = string(M.classes(:)).';
            trVals = M.thresholds(:).';
            % copy values by class name
            for k = 1:numel(keys)
                hit = find(strcmpi(keys{k}, cellstr(trKeys)), 1);
                if ~isempty(hit)
                    vals{k} = trVals(hit);
                end
            end
        end
    catch
        % keep defaults
    end

    thrMap = containers.Map(keys, vals);
end
function tf = maskHas(m, s)
    tf = (s <= size(m,1)) && any(m(s,:));
end
function lab = majorityLabelFromMasks(driftM, delamM, crackM, impactM)
% precedence: Impact > Crack > Delamination > Drift > Normal
    if any(impactM), lab = "Impact";       return; end
    if any(crackM),  lab = "Crack";        return; end
    if any(delamM),  lab = "Delamination"; return; end
    if any(driftM),  lab = "Drift";        return; end
    lab = "Normal";
end
function y = overrideWithDetectors(y, driftM, delamM, crackM, impactM)
% Precedence: Impact > Crack > Delamination > Drift
    S = numel(y);
    for s = 1:S
        if any(impactM(s,:))
            y(s) = "Impact";
        elseif any(crackM(s,:))
            y(s) = "Crack";
        elseif any(delamM(s,:))
            y(s) = "Delamination";
        elseif any(driftM(s,:))
            y(s) = "Drift";
        else
            y(s) = "Normal";
        end
    end
end
function plotConsensusWithClassMarkers(ax, t, Y, driftM, delamM, crackM, impactM)
% Plots mean strain and class-coloured markers with class precedence
% (Impact > Crack > Delamination > Drift) aggregated across sensors.
%
% Colours/markers:
%   Drift:        'o'  [0.85 0.33 0.10]  (orange/red)
%   Delamination: 's'  [0.20 0.62 0.20]  (green)
%   Crack:        'd'  [0.00 0.00 0.00]  (black)
%   Impact:       'p'  [0.49 0.18 0.56]  (purple)

    meanStrain = mean(Y,1);

    % Aggregate "any sensor" masks by class with precedence
    mkI = any(impactM,1);
    mkC = ~mkI & any(crackM,1);
    mkD = ~mkI & ~mkC & any(delamM,1);
    mkR = ~mkI & ~mkC & ~mkD & any(driftM,1);

    % Plot
    cla(ax); hold(ax,'on'); grid(ax,'on');
    hMean = plot(ax, t, meanStrain, 'k-', 'LineWidth', 2, 'DisplayName','Mean strain');

    H = gobjects(1,5);  % legend handles
    H(1) = hMean;

    if any(mkR)
        H(end+1) = plot(ax, t(mkR), meanStrain(mkR), 'o', ...
            'MarkerSize',6, 'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','none', ...
            'LineStyle','none', 'DisplayName','Drift');
    end
    if any(mkD)
        H(end+1) = plot(ax, t(mkD), meanStrain(mkD), 's', ...
            'MarkerSize',6, 'MarkerFaceColor',[0.20 0.62 0.20], 'MarkerEdgeColor','none', ...
            'LineStyle','none', 'DisplayName','Delamination');
    end
    if any(mkC)
        H(end+1) = plot(ax, t(mkC), meanStrain(mkC), 'd', ...
            'MarkerSize',6, 'MarkerFaceColor',[0 0 0], 'MarkerEdgeColor','none', ...
            'LineStyle','none', 'DisplayName','Crack');
    end
    if any(mkI)
        H(end+1) = plot(ax, t(mkI), meanStrain(mkI), 'p', ...
            'MarkerSize',7, 'MarkerFaceColor',[0.49 0.18 0.56], 'MarkerEdgeColor','none', ...
            'LineStyle','none', 'DisplayName','Impact');
    end

    title(ax,'Consensus Mean Strain & Faults');
    xlabel(ax,'Time (s)'); ylabel(ax,'Mean Strain');

    % Clean legend (valid handles only)
    H = H(isgraphics(H));
    legend(ax, H, 'Location','best', 'AutoUpdate','off');
    hold(ax,'off');
end
function maybe_log_scores(score, modelClasses, yTrue, yPred, faultType, fig, useML, useKalman)
% Append per-sensor scores to ./exports/scores_master.csv (robustly).
% Safe to call every run; it will create the file and headers if needed.
%
% Required:
%   score         : S x K matrix from predict()
%   modelClasses  : 1 x K string/char array of class names (order = score columns)
%   yTrue         : S x 1 string or cellstr ground truth labels
%
% Optional (pass [] to skip):
%   yPred         : S x 1 predicted labels (string/cellstr)  [unused but logged if present]
%   faultType     : scalar string/char ("Drift", etc.)
%   fig           : figure handle that may carry SNR via getappdata(fig,'SNR_DB')
%   useML         : logical
%   useKalman     : logical

% ---------- guards ----------
if nargin < 4 || isempty(score) || isempty(modelClasses) || isempty(yTrue)
    return; % nothing to log
end
if nargin < 5, faultType = ""; end
if nargin < 6, fig = []; end
if nargin < 7, useML = true; end
if nargin < 8, useKalman = true; end

% normalize types/shapes
if iscell(yTrue), yTrue = string(yTrue); end
if iscell(yPred), yPred = string(yPred); end
if ischar(modelClasses), modelClasses = string(cellstr(modelClasses)); end
if isrow(modelClasses), modelClasses = modelClasses(:)'; end

[S, ~] = size(score);
modelClasses = lower(strtrim(modelClasses));

% SNR, variant
snrDB = NaN;
if ~isempty(fig) && isgraphics(fig)
    v = getappdata(fig,'SNR_DB');
    if ~isempty(v), snrDB = double(v); end
end
variant = "ML="+string(logical(useML)) + ",Kal=" + string(logical(useKalman));

% Build table row block
T = table();
T.fault   = repmat(string(faultType), S, 1);
T.snrDB   = repmat(snrDB, S, 1);
T.variant = repmat(variant, S, 1);
T.true    = yTrue(:);
if ~isempty(yPred)
    T.pred = yPred(:);
end

% Standardized score_* columns (sorted alphabetically for stability)
allStdNames = lower(["Normal","Drift","Delamination","Crack","Impact"]);
% map each standard class to a column (or zeros if classifier lacks it)
Smat = zeros(S, numel(allStdNames));
for j = 1:numel(allStdNames)
    cls = allStdNames(j);
    idx = find(modelClasses == cls, 1);
    if ~isempty(idx)
        Smat(:,j) = score(:,idx);
    else
        Smat(:,j) = 0; % classifier didn't have this class column
    end
end
for j = 1:numel(allStdNames)
    T.( "score_" + allStdNames(j) ) = Smat(:,j);
end

% Write/append
outDir = fullfile(pwd,'exports');
if ~exist(outDir,'dir'), mkdir(outDir); end
fout = fullfile(outDir, 'scores_master.csv');
if isfile(fout)
    writetable(T, fout, 'WriteMode','append');
else
    writetable(T, fout);
end
end