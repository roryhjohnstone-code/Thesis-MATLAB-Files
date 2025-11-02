function mdl = train_and_cache_model(samplesPerClass, method)
% TRAIN_AND_CACHE_MODEL  Train a classifier and cache it in a stable schema.
% Usage:
%   mdl = train_and_cache_model();                 % default: 300 per class, 'ensemble'
%   mdl = train_and_cache_model(500,'ensemble');   % custom
%
% Produces a struct with fields:
%   .bestClassifier  : classifier object (fitcensemble/fitctree/etc.)
%   .classifier      : alias of bestClassifier
%   .classNames      : string row vector of class labels
%   .featureNames    : string row vector of predictor names
%   .trainedAt       : datetime (UTC)
%   .method          : training method string
%
% Saved to disk (for consistent loads): bestClassifier.mat  (variable: mdl)
% Also cached in memory (appdata): key 'FBG_SHM_ML_MODEL'

    if nargin < 1 || isempty(samplesPerClass), samplesPerClass = 300; end
    if nargin < 2 || isempty(method),          method = 'ensemble';  end

    % 1) Generate training data (features table + categorical labels)
    [trainX, trainY] = generateTrainingData(samplesPerClass);  % table + categorical

    % 2) Train
    mdlStruct = trainClassifier(trainX, trainY, method);

    % 3) Extract the classifier object regardless of wrapper format
    if isfield(mdlStruct,'Model') && ~isempty(mdlStruct.Model)
        clf = mdlStruct.Model;
    elseif isfield(mdlStruct,'bestClassifier') && ~isempty(mdlStruct.bestClassifier)
        clf = mdlStruct.bestClassifier;
    else
        clf = mdlStruct;  % some trainers return the model directly
    end

    % 4) Package a stable schema that runSimulation/ensureClassifier can rely on
    mdl = struct();
    mdl.bestClassifier = clf;
    mdl.classifier     = clf;                                 % legacy alias
    mdl.classNames     = string(categories(trainY))';          % keep real order
    if istable(trainX)
        mdl.featureNames = string(trainX.Properties.VariableNames);
    else
        mdl.featureNames = compose("f%d", size(trainX,2));
    end
    mdl.trainedAt      = datetime('now','TimeZone','UTC');
    mdl.method         = string(method);

    % 5) Cache + save (consistent variable name 'mdl')
    setappdata(0,'FBG_SHM_ML_MODEL', mdl);
    save('bestClassifier.mat','mdl');

    disp('✅ Trained model cached and saved to bestClassifier.mat (variable: mdl)');
end
