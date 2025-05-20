function [beta, intercept, W, T, P, c] = continuumRegression(X, y, alpha, nComponents)
%CONTINUUMREGRESSION Performs Continuum Regression.
%
%   Inputs:
%       X           - Predictor matrix (nSamples x nPredictors)
%       y           - Response vector (nSamples x 1)
%       alpha       - Continuum parameter (0 <= alpha <= 1)
%       nComponents - Number of components to extract
%
%   Outputs:
%       beta        - Regression coefficients vector (nPredictors x 1)
%       intercept   - Intercept term (scalar)
%       W           - Weight matrix (nPredictors x nComponents)
%       T           - Score matrix (nSamples x nComponents)
%       P           - Loading matrix (nPredictors x nComponents)
%       c           - Regression coefficients for components (nComponents x 1)
%
%   Reference:
%       Xie Z, Chen X, Huang G. Optimizing a vector of shrinkage factors 
%       for continuum regression[J]. Chemometrics and Intelligent Laboratory 
%       Systems, 2020, 206: 104141.

% Ensure that X and y are double precision
X = double(X);
y = double(y);

% Step 0: Center the data
X_mean = mean(X, 1);
y_mean = mean(y);
Xc = X - X_mean;
yc = y - y_mean;

[nSamples, nPredictors] = size(Xc);

% Initialize matrices
W = zeros(nPredictors, nComponents);  % Weights
T = zeros(nSamples, nComponents);     % Scores
P = zeros(nPredictors, nComponents);  % Loadings
c = zeros(nComponents, 1);            % Coefficients
Xk = Xc;  % Deflated X
yk = yc;  % Deflated y

for k = 1:nComponents
    % Step 1: Compute the weight vector w_k
    SXX = Xk' * Xk / (nSamples-1);  % Covariance matrix of X
    SXY = Xk' * yk / (nSamples-1);  % Covariance between X and y
    
    % Compute SXX^(- (1 - alpha))
    % Eigenvalue decomposition to compute matrix power
    [eigVectors, eigValues] = eig(SXX);
    eigValues = diag(eigValues);
        
    % Regularization to handle near-zero eigenvalues
    epsilon = 1e-6;
    eigValues(eigValues < epsilon) = epsilon;

    % Compute SXX^(alpha/(1 - alpha)-1)
    SXX_inv_alpha = eigVectors * diag(eigValues .^(alpha/(1-alpha)-1)) * eigVectors';
    
    % Compute weight vector w_k
    wk = SXX_inv_alpha * SXY;
    wk = wk / norm(wk);  % 
    W(:, k) = wk;
    
    % Step 2: Compute the score vector t_k
    tk = Xk * wk;
    T(:, k) = tk;
    
    % Step 3: Compute the loading vector p_k and coefficient c_k
    pk = (Xk' * tk) / (tk' * tk);
    P(:, k) = pk;
    ck = (tk' * yk) / (tk' * tk);
    c(k) = ck;
    
    % Step 4: Deflate X and y with alpha
    Xk = Xk - tk * wk' ;
    yk = yk - ck * tk  ;
end

% Compute regression coefficients beta
beta = W * c;

% Adjust for centering to compute the intercept
intercept = y_mean - X_mean * beta;

end
