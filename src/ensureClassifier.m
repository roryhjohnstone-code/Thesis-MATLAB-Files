function M = ensureClassifier()
% Load trained model saved by ml_train_classifier.m as classifier_ML.mat
% Expects variable M with fields: mdl, mu, sigma, classes.
    if exist('classifier_ML.mat','file') ~= 2
        error('ensureClassifier:NoFile','classifier_ML.mat not found. Train and save the model first.');
    end
    S = load('classifier_ML.mat');   % should contain M
    if ~isfield(S,'M')
        error('ensureClassifier:BadFile','classifier_ML.mat does not contain variable M.');
    end
    M = S.M;
    % light validation
    req = {'mdl','mu','sigma','classes'};
    for k = 1:numel(req)
        if ~isfield(M, req{k})
            error('ensureClassifier:BadStruct','Trained model M missing field: %s', req{k});
        end
    end
end

