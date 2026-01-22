function [outsig] = PHAINmain_TF(insig, corrupted_sgram, mask, param, paramsolver,oracle)

% param
%   .G ........ STFT
%   .G_adj .... inverse STFT
%   .G_diff ... STFT with derivative of window g
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

% define proximal operators
param.proj = @(x) x.*(1-mask) + corrupted_sgram.*mask;
param.prox = @(z) param.thresholding(z);

% first calculate omega from the input signal
omega_y = param.omega(insig);

x_old = insig;
if param.type == "U-PHAIN-TF-GCPA" 

    % setup operator R and its adjoint R*
    phaseCor = exp(-1i*param.phaseCor(omega_y));
    invPhaseCor = exp(1i*param.phaseCor(omega_y));
    param.R = @(z) phaseCor.*z;
    param.R_adj =  @(z) invPhaseCor.*z;

    % combine operators into one: L=DRG, L*=G*R*D*
    param.L = @(x) param.D(param.R(param.G(x)));
    param.L_adj = @(u) param.G_adj(param.R_adj(param.D_adj(u)));
    
    % set starting x, y, z for GCPA
    paramsolver.x0 = insig;
    paramsolver.y0 = zeros(size(param.G(insig)));
    paramsolver.z0 = zeros(size(param.L(insig)));
    
    
    for j = 1:paramsolver.J
        % calculate GCPA
        x_hat = GCPA(param, paramsolver, oracle);
    
        % stopping criterion
        if norm(x_old - x_hat) < paramsolver.epsilon
            break 
        end
    
        % calculate new IF and redefine L and L_adj
        omega_x_hat = param.omega(x_hat);
    
       % redefine R, R*, L and L* with new instFreq
        phaseCor = exp(-1i*param.phaseCor(omega_x_hat));
        invPhaseCor = exp(1i*param.phaseCor(omega_x_hat));
        param.R = @(z) phaseCor.*z;
        param.R_adj =  @(z) invPhaseCor.*z;

        param.L = @(x) param.D(param.R(param.G(x)));
        param.L_adj = @(u) param.G_adj(param.R_adj(param.D_adj(u)));
        
    
        x_old = x_hat;
    end
    % output spectrogram
    outsig = param.proj(param.G(x_hat));

end
