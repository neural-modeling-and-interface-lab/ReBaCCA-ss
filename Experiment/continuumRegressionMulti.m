function [A, B, U, V, P, Q,objValueTotal,objComponent,objCompTrue,varExplained,varExplainedX,varExplainedY] =...
    continuumRegressionMulti(X, Y, alpha, nComponents, tol, maxIter,percentVar)
%CONTINUUMREGRESSIONMULTIY Performs Continuum Regression with multi-dimensional Y.

%   Inputs:
%       X           - Predictor matrix (nSamples x nPredictors)
%       Y           - Response matrix (nSamples x nResponses)
%       alpha       - Continuum parameter (0 <= alpha <= 1)
%       nComponents - Number of components to extract
%       tol         - Tolerance for convergence (e.g., 1e-6)
%       maxIter     - Maximum number of iterations for convergence (e.g., 500)
%
%   Outputs:
%       A           - Weight matrix for X (nPredictors x nComponents)
%       B           - Weight matrix for Y (nResponses x nComponents)
%       U           - Score matrix for X (nSamples x nComponents)
%       V           - Score matrix for Y (nSamples x nComponents)
%       P           - Loading matrix for X (nPredictors x nComponents)
%       Q           - Loading matrix for Y (nResponses x nComponents)
%
%   Reference:
%       Extension of Continuum Regression to multi-dimensional Y.

% Ensure that X and Y are double precision
X = double(X);
Y = double(Y);

% Step 0: Center and normalize the data
X_mean = mean(X, 1);
Y_mean = mean(Y, 1);
Xc = X - X_mean;
Yc = Y - Y_mean;


% a=max(eig(cov(Xc)));
% b=max(eig(cov(Yc)));
% Xc=Xc/sqrt(a);
% Yc=Yc/sqrt(b);


[nSamples, nPredictors] = size(Xc);
[~, nResponses] = size(Yc);

% Initialize matrices
A = zeros(nPredictors, nComponents);  % Weights for X
B = zeros(nResponses, nComponents);   % Weights for Y
U = zeros(nSamples, nComponents);     % Scores for X
V = zeros(nSamples, nComponents);     % Scores for Y
P = zeros(nPredictors, nComponents);  % Loadings for X
Q = zeros(nResponses, nComponents);   % Loadings for Y
Xk = Xc;  % Deflated X
Yk = Yc;  % Deflated Y

% Initialize v_k according to alpha value
if alpha < 0.5
    [Acca,Bcca] = canoncorr(Xk,Yk);
    wk_initial=Acca;
    vk_initial=Bcca;
else 
    [Px,~]=pca(Xk);
    [Py,~]=pca(Yk);
    wk_initial=Px;
    vk_initial=Py;
end

objValueTotal=[];
objComponent=zeros(3,nComponents);
objCompTrue=zeros(3,nComponents);
varExplainedX=zeros(1,nComponents);
varExplainedY=zeros(1,nComponents);
stoppedIdx=nComponents;

totalVarianceX = sum(diag(cov(X)));
totalVarianceY = sum(diag(cov(Y)));

for k = 1:nComponents

    objValue=0;

    %% Initialize v_k 
    vk = vk_initial(:,k);
    % vk = randn(nResponses,1);
    vk = vk / norm(vk);
    targetYk=Yk*vk;
    
    wk=wk_initial(:,k);
    % wk= randn(nPredictors,1);
    wk=wk / norm(wk);

    % Iterative algorithm to find w_k and v_k
    for iter = 1:maxIter
        % Store previous w_k and v_k for convergence check
        vk_old = vk;
        wk_old = wk;

        [~,~,wk]=continuumRegression(Xk,targetYk,alpha,1);  
        targetXk=Xk*wk;

        [~,~,vk]=continuumRegression(Yk,targetXk,alpha,1);
        targetYk=Yk*vk;

        objValue=[objValue (wk'*Xk'*Xk*wk/(wk'*wk)/(nSamples-1)/totalVarianceX)^alpha*...
            ((wk'*Xk'*Yk*vk)^2/(wk'*Xk'*Xk*wk)/(vk'*Yk'*Yk*vk))^(1-alpha)*...
            (vk'*Yk'*Yk*vk/(vk'*vk)/(nSamples-1)/totalVarianceY)^alpha ];
        
        % In terms of the original data
        % objValue=[objValue (wk'*Xc'*Xc*wk/(wk'*wk)/(nSamples-1))^alpha*...
        %     ((wk'*Xc'*Yc*vk)^2/(wk'*Xc'*Xc*wk)/(vk'*Yc'*Yc*vk))^(1-alpha)*...
        %     (vk'*Yc'*Yc*vk/(vk'*vk)/(nSamples-1))^alpha ];

        % Check for convergence
        % if norm(vk - vk_old)/norm(vk_old) < tol && norm(wk - wk_old)/norm(wk_old) < tol
        if (objValue(end)-objValue(end-1))/objValue(end-1) < tol
            vk = vk_old;
            wk = wk_old;
            objValueTotal=[objValueTotal objValue(end-1)];
            break;
        end
    end
    
    if iter == maxIter
        warning('Maximum iterations reached without convergence at component %d.', k);
        objValueTotal=[objValueTotal objValue(end)];
    end
            
    objComponent(1,k)=(wk'*Xk'*Xk*wk/(wk'*wk)/(nSamples-1)) / totalVarianceX;
    objComponent(2,k)=((wk'*Xk'*Yk*vk)^2/(wk'*Xk'*Xk*wk)/(vk'*Yk'*Yk*vk));
    objComponent(3,k)=(vk'*Yk'*Yk*vk/(vk'*vk)/(nSamples-1)) / totalVarianceY;
    
    objCompTrue(1,k)=objComponent(1,k)^alpha;
    objCompTrue(2,k)=objComponent(2,k)^(1-alpha);
    objCompTrue(3,k)=objComponent(3,k)^alpha;
    
    % % 2025-03-25 Normalize based on alpha
    % scale = 1 / ((1-alpha) * std(Xk * wk) + alpha);
    % wk = wk * scale;
    % scale = 1 / ((1-alpha) * std(Yk * vk) + alpha);
    % vk = vk * scale;

    A(:, k) = wk;
    B(:, k) = vk;
    
    % Step 3: Compute scores t_k and u_k
    tk = Xk * wk;
    uk = Yk * vk;
    U(:, k) = tk;
    V(:, k) = uk;
    
    % Step 4: Compute loadings p_k and q_k
    pk = (Xk' * tk) / (tk' * tk);
    P(:, k) = pk;
    qk = (Yk' * uk) / (uk' * uk);
    Q(:, k) = qk;
    
    varExplainedX(k)=objComponent(1,k);
    varExplainedY(k)=objComponent(3,k);

    % Step 5: Deflate X and Y
    Xk = Xk - tk * pk';
    Yk = Yk - uk * qk';

    % 2025-03-20 Make deflation with alpha
    % Xk = Xk - tk * wk' / std(tk)^(1-alpha);
    % Yk = Yk - uk * vk' / std(uk)^(1-alpha);
    
    if sum(geomean([varExplainedX; varExplainedY]))>=percentVar
        stoppedIdx=k;
        break;
    end
    
    % Xk = Xk - tk * wk';
    % Yk = Yk - uk * vk';
end

varExplained= geomean([varExplainedX; varExplainedY]);

A=A(:,1:stoppedIdx);
B=B(:,1:stoppedIdx);
U=U(:,1:stoppedIdx);
V=V(:,1:stoppedIdx);
P=P(:,1:stoppedIdx);
Q=Q(:,1:stoppedIdx);
objComponent=objComponent(:,1:stoppedIdx);
objCompTrue=objCompTrue(:,1:stoppedIdx);
varExplained=varExplained(:,1:stoppedIdx);
varExplainedX=varExplainedX(:,1:stoppedIdx);
varExplainedY=varExplainedY(:,1:stoppedIdx);

end
