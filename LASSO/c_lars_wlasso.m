function [Beta,Asupp,Lamda] = c_lars_wlasso(y,X,K,w,intercept,scaleX)
% This algorithm refers to the paper of Muhammad Naveed Tabassum and Esa Ollila
% "Sequential adaptive elastic net approach for single-snapshot source
% localization" J. Acoust. Soc. Am. / 22 May 2018

% Inputs:
% y: the n complex-valued observations 
% X: the design matrix
% w: the weight vector
% K: the level of sparsity 
% intercept: the indicator of if get y centerted or not
% scaleX: the indicator of if get X normalized or not 

% Outputs:
% Beta: the beta coefficients at non-zero knots
% Asupp: the index collection of non-zero predictors
% Lamda: the konts where the active set changes

K=K+1;
[n,p] = size(X);
if nargin < 6
    scaleX = 1; 
end
if nargin < 5
    intercept = 1; 
end
if nargin < 4
    w = ones(p,1);
end
if nargin < 3
    K = min(n,p);
end

% Get X normailized
if scaleX
    % Default: scale X to have column norms eq to sqrt(n) (n = nr of rows)
    sdX = sqrt(sum(X.*conj(X)));
    X = bsxfun(@rdivide, X, sdX);
end
% Get y centered
if intercept
    meanX = mean(X); % Calculate the mean of each column of X,return a vector
    meany= mean(y);
    X = bsxfun(@minus, X, meanX); % bsxfun(): Apply element-wise operation to two arrays
    y = y-meany;
end


%% initialize 
beta = zeros(p,1);
Beta = zeros(p,K);
Lamda=zeros(K,1);
X = X./repmat(w',n,1); % X = X*inverse(diag(w))
r = y-X*beta; % Calculate the first residual vector
cor = X'*r;
[lamda,Asupp] = max(abs(cor)); % Record the first active predictor and corresbonding index
Lamda(1)=lamda;
Delta = zeros(p,1);

%% Loop for calculating betas and corresbonding lamdas
% We have a non-zero coefficient already, so we start from k=2:K
for k=2:K
    delta = (1/Lamda(k-1))*(X(:,Asupp)\r);  % "\" is the inverse operator of matlab
    Delta(Asupp)=delta; % Record the direction of lars
    notA=setdiff(1:p,Asupp); % Find the index of non-active predictors
    notact=length(notA);
    Gamma=zeros(1,p); % Record gammas of non-active predictors
    for i = 1:notact
        notindex = notA(i);
        % We find the knot lamda(k) by solving all unary quadratic equation
        % of non-active predictors
        cl = X(:,notindex)'*r;
        bl = X(:,notindex)'*(X(:,Asupp)*delta);
        Ax2 = abs(bl)^2-1; % the coeffcient of quadratic item
        Bx1 = 2*Lamda(k-1)-2*real(cl*conj(bl)); % the coefficent of first-order item
        Cx0 = abs(cl)^2-Lamda(k-1)^2; % the coefficient of the constant item
        uqp = [Ax2, Bx1, Cx0]; % Build the Unary quadratic equation
        gammas=roots(uqp); % Find the roots
        
        if all(gammas > 0) % The all function reduces such a vector of logical values to a single condition
            gamma = min(gammas);
        else
            gamma = subplus(max(gammas));
        end
        
        if ~isreal(gammas) % if these roots contain complex values, set gamma equal to zero
            gamma = 0;
        end
        
        Gamma(notindex)=gamma; % Record the stepsize of the direction of active predictors
        cor(notindex)=cl-gamma*bl; % the correlation of the next active predictor is changed
    end
    [mingamma,pos] = min(Gamma(notA)); % Find the smallest gamma 
    actindex = notA(pos); %  Find the index of the smallest gamma 
    Lamda(k) = Lamda(k-1) - mingamma; % Update the penalty coefficient
    Asupp =[Asupp,actindex]; % Add the new predictor to the active set
    beta = beta+mingamma*Delta; % Update the coefficient vector
    Beta(:,k)=beta; % Record the coefficient vector
    r = y-X*beta; % Update the residual
end

%% Solve the original coefficients
Asupp = Asupp(1:K-1)';
Lamda = Lamda(1:K-1)';
Beta = Beta./repmat(w,1,K); % Divide the weights 
if scaleX % Compensate the scaling and centering
    Beta = bsxfun(@rdivide,Beta,sdX');
end
if intercept
    Beta = [meany - meanX*Beta;  Beta];
end
Beta=Beta(:,2:K);
    
    
    
    
    
    