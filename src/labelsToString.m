function s = labelsToString(x)
% Robustly convert class labels (ClassLabel, categorical, cellstr, char) to string row vector.
    if isa(x, 'classreg.learning.internal.ClassLabel')
        % Convert each ClassLabel to char, then to string
        s = string(arrayfun(@char, x, 'UniformOutput', false));
    elseif iscategorical(x)
        s = string(x);
    elseif iscellstr(x)
        s = string(x);
    elseif isstring(x)
        s = x;
    elseif ischar(x)
        s = string(cellstr(x)); % 1×N char row -> 1×1 string
    else
        % Last-resort: try string() on a cell array of chars
        try
            s = string(x);
        catch
            s = string(arrayfun(@char, x, 'UniformOutput', false));
        end
    end
    s = s(:).';  % make it a row vector
end
