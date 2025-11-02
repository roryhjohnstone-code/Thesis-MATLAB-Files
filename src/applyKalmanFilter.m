function y = applyKalmanFilter(x, q, r)
% applyKalmanFilter  Minimal 1D scalar Kalman smoother for a single channel.
% y = applyKalmanFilter(x, q, r)
%   x : row or column vector (signal)
%   q : process noise variance   (e.g., 0.05)
%   r : measurement noise var.   (e.g., 1e-4)

x = x(:);                  % column
n = numel(x);
y = zeros(n,1);

% simple constant-state model
A = 1; H = 1;
Q = q; R = r;
xhat = x(1);               % initial state
P    = 1;                  % initial covariance

for k = 1:n
    % Predict
    xhat = A * xhat;
    P    = A * P * A' + Q;

    % Update
    K    = P * H' / (H * P * H' + R);
    xhat = xhat + K * (x(k) - H * xhat);
    P    = (1 - K * H) * P;

    y(k) = xhat;
end

% return in original shape (row/col)
if isrow(x.'), y = y.'; end
end