function populateFaultTable(fig, time, driftDetected, ~, delaminationDetected, crackDetected, impactDetected)
% Populate the Fault Summary uitable (one row per sensor, first event only).
% Signature kept backward-compatible with your existing calls.

    if nargin < 6, crackDetected = false(size(driftDetected)); end
    if nargin < 7, impactDetected = false(size(driftDetected)); end
    if isempty(fig) || ~isgraphics(fig), return; end

    % Find or create the summary table (top-right)
    ht = getappdata(fig,'tblFaultSummary');
    if isempty(ht) || ~isvalid(ht)
        % Create a reasonable default position (matches your GUI layout)
        ht = uitable('Parent',fig,'Units','normalized',...
            'Position',[0.73 0.63 0.24 0.27], ...
            'ColumnName',{'Sensor','Fault Type','Start Time (s)','End Time (s)'});
        setappdata(fig,'tblFaultSummary',ht);
    else
        set(ht,'ColumnName',{'Sensor','Fault Type','Start Time (s)','End Time (s)'});
    end

    % Ensure masks are S×N
    [S,N] = size(driftDetected); %#ok<ASGLU>
    if size(driftDetected,2) ~= numel(time), driftDetected = driftDetected.'; end
    if size(delaminationDetected,2) ~= numel(time), delaminationDetected = delaminationDetected.'; end
    if size(crackDetected,2) ~= numel(time), crackDetected = crackDetected.'; end
    if size(impactDetected,2) ~= numel(time), impactDetected = impactDetected.'; end
    S = size(driftDetected,1);

    % Precedence per sensor: Impact > Crack > Delam > Drift
    rows = cell(0,4);
    for s = 1:S
        mImp = logical(impactDetected(s,:));
        mCrk = logical(crackDetected(s,:));
        mDel = logical(delaminationDetected(s,:));
        mDr  = logical(driftDetected(s,:));

        [lab,t1,t2] = pickFirst(time, mImp, mCrk, mDel, mDr);
        if ~isempty(lab)
            rows(end+1,:) = {sprintf('Sensor %d',s), char(lab), t1, t2}; %#ok<AGROW>
        end
    end

    if isempty(rows)
        set(ht,'Data',cell(0,4));
    else
        % uitable wants char/double, not string
        for r = 1:size(rows,1)
            if ~ischar(rows{r,1}), rows{r,1} = char(rows{r,1}); end
            if ~ischar(rows{r,2}), rows{r,2} = char(rows{r,2}); end
            rows{r,3} = double(rows{r,3});
            rows{r,4} = double(rows{r,4});
        end
        set(ht,'Data',rows);
    end
end

% ---- helpers ----
function [lab,t1,t2] = pickFirst(t, mImp, mCrk, mDel, mDr)
    lab=''; t1=[]; t2=[];
    if any(mImp), [t1,t2] = firstSpan(t,mImp); lab='Impact';       return; end
    if any(mCrk), [t1,t2] = firstSpan(t,mCrk); lab='Crack';        return; end
    if any(mDel), [t1,t2] = firstSpan(t,mDel); lab='Delamination'; return; end
    if any(mDr),  [t1,t2] = firstSpan(t,mDr);  lab='Drift';        return; end
end

function [t1,t2] = firstSpan(t, m)
    m = logical(m(:).'); d = diff([0 m 0]);
    b = find(d==1,1,'first'); e = find(d==-1,1,'first')-1;
    if isempty(b), t1=[]; t2=[]; else, t1=t(b); t2=t(e); end
end

