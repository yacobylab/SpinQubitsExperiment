function filt = updateKalman(filt,grad,gradDev)
Q=diag([5 1 1]); % Process noise, half-assed guess.
H = [1 0 0]; % Perform a correction step.
% P: the P matrix 
% F: the update matrix 
% We bring in: x, F, P 
% We update: P, x

x_priori = filt.F * filt.x ; % Predict the new state
P_priori = filt.F * filt.P * (filt.F') + Q; % Predict new process matrix.
grad = abs(grad)*sign(x_priori(1));  % Guess sign of gradient.
filt.y = grad - x_priori(1); % Innovation
S = P_priori(1,1) + gradDev^2; % Innovation covariance
K = P_priori * (H') / S; % Optimal gain
filt.x = x_priori + K * filt.y; % a-posteriori value of x
filt.P = (eye(size(filt.P)) - K * H) * P_priori; % New process matrix
end