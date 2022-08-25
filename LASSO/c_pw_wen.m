function [Beta_coeff,Best_alpha,Active_set,min_RSS] = c_pw_wen(y,X,w,K,Alpha,DEBIAS)
% This algorithm refers to the paper of Muhammad Naveed Tabassum and Esa Ollila
% "Sequential adaptive elastic net approach for single-snapshot source
% localization" J. Acoust. Soc. Am. / 22 May 2018

% Inputs:
% y: the n complex-valued observations 
% X: the design matrix
% w: the weight vector
% K: the level of sparsity 
% Alpha: the dense grid of 0<alpha<=1
% DEBIAS: Get beta debiased,only keep the non-zero elements

% Outputs:
% Beta_coeff: the sparse coefficient vector with sparsity K
% Best_alpha: the alpha with minimum RSS
% Active_set: the index collection of non-zero predictors
% min_RSS: the smallest RSS

[n,p] = size(X);
if nargin < 6
    DEBIAS = 1;
end
if nargin < 5
    Alpha = 1;
end
if nargin < 4
    K=min(n,p);
end
if nargin < 3
    w=ones(p,1);
end

%% Initialize
nral = length(Alpha);
[Beta,Asupp,Lamda] = c_lars_wlasso(y,X,K,w,0,0); % Solve the initial values of Beta, Asupp and Lamda
mLamda = zeros(K,nral); % Store penalty items
mLamda(:,1) = Lamda;
mAsupp = zeros(K,nral); % Store active sets
mAsupp(:,1) = Asupp;
mBeta = zeros(p, nral); % Store beta coefficients
mBeta(:,1) = Beta(:,K); 
ya = [y',zeros(1,p)]'; % The augmented form of y
mBetaLS = zeros(K,nral); % Store the unbiased beta coefficients
mBetaLS(:,1) = X(:,Asupp)\y; % Calculate the first unbiased beta coefficients
mRSS = zeros(nral,1); % Store residuals
r1 = y-X(:,Asupp)*mBetaLS(:,1);
mRSS(1) = sum(r1'*r1);

%% Start Loops to find out beta,alpha and Asupp 
for i = 2:nral
    % This is the loop to search the dense grid of alpha
    nK = mLamda(K,i-1)*(1-Alpha(i));
    Xa = [X',sqrt(nK)*eye(p)]'; % The augmented form of X
    [BetaK,AsuppK,pKnots] = c_lars_wlasso(ya,Xa,K,w,0,0); % Solve the K non-zero solution at specific alpha
    Knots = pKnots/Alpha(i); % Calculate the Kth knot of lamda at specific alpha
    mLamda(:,i) = Knots;
    mAsupp(:,i) = AsuppK;
    mBeta(:,i) = BetaK(:,K);
    BetaLS = X(:,AsuppK)\y;
    mBetaLS(:,i) = BetaLS;
    ri = y-X(:,AsuppK)*BetaLS;
    RSSi = sum(ri'*ri);
    mRSS(i) = RSSi;
end
    
%% Find out the beta solution with samllest RSS
[min_RSS,Alindex] = min(mRSS); % Find the alpha with minimum residual
Best_alpha = Alpha(Alindex);
Active_set = mAsupp(:,Alindex);
Beta_coeff = mBeta(:,Alindex);
Beta_coeff = (ones(p,1)+mLamda(K,Alindex)*(1-Best_alpha)*w).*Beta_coeff; % Correction for double shrinkage

if DEBIAS
    Beta_coeff = mBetaLS(:,Alindex); % Get beta unbiased
end
    





