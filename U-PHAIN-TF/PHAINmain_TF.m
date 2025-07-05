function [outsig] = PHAINmain_TF(insig, gapped_spec, mask, param, paramsolver,oracle)

% param
%   .S ........ STFT
%   .S_adj .... inverse STFT
%   .S_diff ... STFT with derivative of window g
%   .omega .... instantaneous frequency
%   .R ........ phase correction
%   .R_adj .... adjoint of phase correction
%   .D ........ time-directional variation
%   .D_adj .... adjoint of time-directional variation
%   .type ..... type of algorithm


% paramsolver
%   .sigma .... step size
%   .tau ...... step size
%   .eta ...... step size (only for GCPA)
%   .alpha .... relaxation paramter
%   .lambda ... threshold
%   .epsilon .. stop criterion
%   .x0 ....... initial value of primal variable
%   .y0 ....... initial value of dual variable
%   .z0 ....... initial value of second dual variable (only for GCPA)
%   .I ........ number of inner iterations
%   .J ........ number of outer iterations



%%

% define proximal functions
soft = @(z, lambda) sign(z).*max(abs(z) - lambda, 0);
param.proj = @(x) x.*(1-mask) + gapped_spec.*mask;
param.prox = @(z) soft(z, paramsolver.lambda/paramsolver.sigma);

% first calculate omega from the input signal
omega_y = param.omega(insig);

x_old = insig;
if param.type == "U-PHAIN-TF-GCPA" 
    % define the analysis and synthesis operators
    param.L = @(x) param.D(param.R(param.S(x), omega_y));
    param.L_adj = @(u) param.S_adj(param.R_adj(param.D_adj(u), omega_y));
    
    % set starting x, y, z for GCPA
    paramsolver.x0 = insig;
    paramsolver.y0 = zeros(size(param.S(insig)));
    paramsolver.z0 = zeros(size(param.L(insig)));
    
    
    for j = 1:paramsolver.J
        % calculate GCPA
        x_hat = GCPA(param, paramsolver,oracle);
    
        % stopping criterion
        if norm(x_old - x_hat) < paramsolver.epsilon
            break 
        end
    
        % calculate new IF and redefine L and L_adj
        omega_x_hat = param.omega(x_hat);
    
        param.L = @(x) param.D(param.R(param.S(x), omega_x_hat));
        param.L_adj = @(u) param.S_adj(param.R_adj(param.D_adj(u), omega_x_hat));
    
        x_old = x_hat;
    end
    % output spectrogram
    outsig = param.proj(param.S(x_hat));

end
