function [Beta_coeff,Active_set] = SeqAEN(y,X,K,Alpha,Mul)
% This algorithm refers to the papers of Muhammad Naveed Tabassum and Esa Ollila
% "Sequential adaptive elastic net approach for single-snapshot source
% localization" J. Acoust. Soc. Am. / 22 May 2018 AND
% "Single-snapshot DoA Estimation using Adaptive Elastic Net in the Complex Domain" 
% 4th Int. Workshop on Compressed Sensing on Radar, Sonar, and Remote Sensing (CoSeRa 2016)

% Inputs:
% y: the n complex-valued observations 
% X: the design matrix
% Alpha: the dense grid of 0<alpha<=1
% DEBIAS: Get beta debiased,only keep the non-zero elements
% Mul: the size of initial solution

% Outputs:
% Beta_coeff: the final coefficient vector with sparsity order K
% Active set: the collection of the indexes of non-zero predictors
% Lamda: the collection of lamdas at 3K, 2K and K

[n,p] = size(X);
if nargin < 5
    Mul = 3;
end
if nargin < 4
    Alpha = 1;
end
if nargin < 3
    K = min(n,p);
end

if Mul*K > min(n,p)
    fprintf('K is too large, it should be smaller than min(nrowX,ncolX)/Mul\n')
    K = floor(min(n,p)/Mul);
end

%% Initialize
Alvec = zeros(Mul,1); % Record the alpha with minimum RSS at 3 different sparsity order
RSSvec = zeros(Mul,1); % Record the minimum RSS at 3 different sparsity order
w = ones(p,1); % Let the initial weight equal to 1
X_iK = X; % Record the active predictors in each loop
w_iK = w; % Reacord the weights in eahc loop
Asupp_iK = 1:p;  % Record the indexes of active predictors in each loop

%% Perform c-pw-wen for three times to reach the final solution with sparsity order K
for i = Mul:-1:1
    [Beta_iK,alpha_iK,Active_iK,RSS_iK] = c_pw_wen(y,X_iK,w_iK,i*K,Alpha,i==1);
    Asupp_iK = Asupp_iK(Active_iK);
    Alvec(i) = alpha_iK;
    RSSvec(i) = RSS_iK;
    if i == 1
        Beta_coeff = Beta_iK;
        Active_set = Asupp_iK;
        break
    end
    Beta_iK = Beta_iK(Active_iK);
    X_iK = X_iK(:,Active_iK);
    w_iK = 1./abs(Beta_iK);
end










