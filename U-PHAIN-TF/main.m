%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%         PHASE-AWARE AUDIO INPAINTING in TF domain (U-PHAIN-TF)          %
%      Variant with the generalized Chambolle-Pock algorithm (GCPA)       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear 
clc
ltfatstart
addpath(genpath('../dataset/'))
%% loading

soundDir = "../dataset/DPAI_originals/";
ext = ".wav";
Sounds = dir(soundDir + "*" + ext);
NN = length(Sounds);
data = cell(NN,1);
info = audioinfo(Sounds(1).name);
fs = info.SampleRate;
for nn = 1:NN
    data{nn} = audioread(Sounds(nn).name);
end
clear audio info

%% define solution and tables
methodLabels = {'U_PHAIN_TF_GCPA'};
for i = 1:length(methodLabels)
    solution.(methodLabels{i}) = cell(NN);  % initialization of restored results
end
SNR_sgram  = NaN(NN, 6, length(methodLabels)); % 6 masks, NN=number of sigs
SNR_sig  = NaN(NN, 6, length(methodLabels));
TIME = NaN(NN, 6, length(methodLabels));
REC_sig = cell(NN, 6, length(methodLabels));
%% parameters

% parameter settings for STFT/DGT
w = 2048; % window length
a = w/4; % hop size
M = w; % number of freq. rows
wtype = 'hann'; % window type
phasetype = 0; % 0-freqinv, 1-timeinv

% setup tight window and its derivative 
g = gabtight(wtype, a, M, w);

if wtype == "hann"
    % derivative of Hann window
    x = (0:w-1)'/(w);
    g_diff = -0.5*sin(2*pi.*x)*max(g);
else 
    % for other window functions
    g_diff = numericalDiffWin(g);
end

% setup DGT, invDGT and derivative of DGT using LTFAT
param.G = @(x) comp_sepdgtreal(x, g, a, M, phasetype);
param.G_adj = @(u) comp_isepdgtreal(u, g, size(u,2)*a, a, M, phasetype);
param.G_diff = @(x) comp_sepdgtreal(x, g_diff, a, M, phasetype);

% definition of instantaneous frequency (omega)
param.omega = @(x) calcInstFreq(param.G(x), param.G_diff(x), M, w);

% define phase rotations and time-directional difference (D)
param.phaseCor = @(omega) rotations(omega, a, M);
param.D = @(z) -diff(z,1,2);
param.D_adj = @(z) [z(:,1), (diff(z,1,2)), -z(:,end)];

% minimal padding around each spectrogram gap
pad = 4;

% settings for generalized Chambolle-Pock algorithm
paramsolver.epsilon = 0.001;  % for stopping criterion
paramsolver.tau = 0.25;  % step size
paramsolver.sigma = 1;  % step size
paramsolver.eta = 4; % step size
paramsolver.alpha = 1;  % relaxation parameter
paramsolver.lambda = 0.01; % threshold (regularization parameter)

% p = 0.9; % uncomment for p-shrinkage thresholding
% alpha = 10e-2; % uncomment for smooth-hard threshodling 

% set inner and outer iterations for U−PHAIN
paramsolver.I = 500;
paramsolver.J = 10;

% thresholding operators (default soft thresholding)
soft = @(z, lambda) sign(z).*max(abs(z) - lambda, 0);
% smoothHard = @(z, lambda) (abs(z) >= lambda).*(z.*exp(-alpha./(exp(abs(z)-lambda)-1).^2));
% pshrink = @(z, lambda) sign(z).*max(abs(z) - lambda^(2-p)*abs(z).^(p-1), 0);
% energyNorm = @(z, lambda) z/(1 + lambda/2);
% l2norm = @(z, lambda) max(1 - lambda/norm(z,2), 0)*z;

% set thresholding (for soft */paramsolver.sigma* can be ommited)
param.thresholding = @(z) soft(z, paramsolver.lambda/paramsolver.sigma);

%% get masks from dpai
masks = dir("../spectrogram_masks/*.mat");
num_of_masks = length(masks);
my_masks = zeros(num_of_masks,157);
for i=1:num_of_masks
    my_var = "C"+i;
    tmp = load("../spectrogram_masks/"+masks(i).name).(my_var);
    my_masks(i,:) = tmp;
end

%% data preparation
% reshape data, so its length is divisible by lcm(a,M) (for LTFAT)
div = lcm(a,M);
for i=1:length(data)
    datasize = size(data{i},1);
    remainder = mod(datasize,div);
    if remainder ~= 0
        % depending on the amount of remainder shorten or expand signal
        if remainder > div/2
            padding = zeros(div-remainder,1);
            data{i} = [data{i}; padding];
        else
            data{i} = data{i}(1:datasize-remainder);
        end
    end
end

%% testing

for nn=1:NN % iterate signals

    current_signal = data{nn};
    sgram = param.G(current_signal); % compute sgram from current signal
    
    for m=1:size(my_masks,1) % iterate masks 1--6 from JanssenTF

        % shorten masks if incompatible length with sgram
        if size(sgram,2)~= size(my_masks(m,:),2)
            my_masks = my_masks(:,1:size(sgram,2));
            disp("reduced size of masks")
        end

        current_mask  = my_masks(m,:);
        corrupted_sgram = sgram.*current_mask;

        % save corrupted spectrogram as solution
        solution.(methodLabels{1}){nn,m} = corrupted_sgram;
        
        % find indexes with zeros
        indxs_gaps = find(all(corrupted_sgram == 0,1));
        tic
        for gap=1:5 % iterate gaps
            % get indexes of gap
            gap_indx = indxs_gaps(gap*m-m+1:gap*m);
            
            % the first gap_index-pad-1 must be divisible by w/a
            [l_pad, r_pad] = fix_pad(pad, gap_indx, a, w);
            segment_indx = (gap_indx(1)-l_pad:gap_indx(end)+r_pad);
            
            % in case bigger padding is used
            if segment_indx(end)>size(sgram,2)
                segment_indx = segment_indx(:,1:end-w/a);
            elseif segment_indx(1)<= 0
                segment_indx = segment_indx(:, 1+w/a:end);
            end

            % get cutout spectrogram, mask, signal
            segment.mask = current_mask(segment_indx);
            segment.corrupted = corrupted_sgram(:,segment_indx);
            segment.corrupted_signal = param.G_adj(segment.corrupted);
            

            % normalize input signal and cutout spectrogram
            segment.max = max(abs(segment.corrupted_signal));
            segment.n_oracle = sgram(:,segment_indx)/segment.max;
            segment.n_signal = segment.corrupted_signal/segment.max;
            segment.n_sgram = segment.corrupted/segment.max;
            
            param.type = "U-PHAIN-TF-GCPA";
            % get reconstructed segment
            [segment.solution] = PHAINmain_TF(segment.n_signal, segment.n_sgram, segment.mask, param, paramsolver,segment.n_oracle);
            % save rec. segment to solution
            solution.(methodLabels{1}){nn,m}(:,segment_indx) = segment.solution*segment.max;

        end
        TIME(nn,m,1) = toc;
        SNR_sgram(nn,m,1) = snr(sgram,sgram-solution.(methodLabels{1}){nn,m});

        signal_data = param.G_adj(sgram);
        signal_result = param.G_adj(solution.(methodLabels{1}){nn,m});

        % write audio to file
        % name_audio = "results/example"+nn+"mask"+m+".wav";
        % audiowrite(name_audio,signal_result,fs);

        REC_sig{nn,m,1} = signal_result;
        SNR_sig(nn,m,1) = snr(signal_data,signal_data-signal_result);
        fprintf('done for example: %d, mask: %d!\n',nn,m)
    end
end
%% plot

% note that the results in the paper are computed from .wav signals
masks = 1:6;
SNRs = mean(SNR_sig);
figure;
plot(masks, SNRs);
xlabel("gap length");
ylabel("mean SNR [dB]")

%% functions

%   [1] Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa,
%       "Representation of complex spectrogram via phase conversion,"
%       Acoustical Science and Technology, vol.40, no.3, May 2019. (Open Access)
%   The below functions were adapted from [1]

function IF = calcInstFreq(spec,diffSpec,fftLen,winLen,flooringCoeff)
    % calc_IF: Calculating instantaneous frequency at each bin.
    
    if ~exist('flooringCoeff','var') || isempty(flooringCoeff)
        flooringCoeff = 1e-10;
    end
    
    powSpec = abs(spec).^2; % power spectrogram
    powSpec = powSpec + flooringCoeff*max(powSpec(:)); % avoiding division by zero
    
    IF = -imag(diffSpec.*conj(spec)./powSpec); % calculating IF by Eq. (21) of [1]
    IF = (fftLen/winLen)*IF; % compensation necessary when "fftLen ~= winLen"
end

function phaseCor = rotations(IF,shiftLen,fftLen)
    % instPhaseCorrection: Calculating instantaneous-phase-corrected spectrogram.
    
    sigLen = shiftLen*size(IF,2); % L (= a * N) : signal length
    freqShift = sigLen/fftLen;    % b (= L / M) : frequency stepsize
    
    idxVariation = freqShift*IF*shiftLen/sigLen;   % b * delta * a / L (in Eq. (29) of [1])
    phaseCor = 2*pi*mod(cumsum(idxVariation,2),1); % mod for avoiding huge value
end
