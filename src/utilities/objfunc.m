function [L,dLdb] = objfunc(b,X,y)
    % For given (fixed) predictor values X, compute the probabalility of
    % observing 1, i.e. Pr(y=1), as function of the coefficients b. The log
    % likelihood function to be maximized is:
    %  
    %    L = sum(y_i*log(p(y=1)) + (1-y_i)*log(1-p(y=1)),
    %                       i = 1 ... n
    %
    % where p(y=1) = 1 / (1 + exp(-b'*Xi)
    %
    % Inputs:
    %
    %   b    Kx1 vector of unknown coefficients
    %   X    X = [X1 X2 ... XK] NxK matrix of K predictors, with N observations each
    %   y    Nx1 vector of responses
    
    [N,K] = size(X);
    if numel(b) ~= K
        error('Vector of coefficients, b, must be of size Kx1, where K is the number of predictors.')
    end
    
    if numel(y) ~= N
        error('Vector of responses, y, must be of size Nx1, where N is the number of observations.')
    end
    
    % Probabilities of y = 1 for all observations
    p = 1 ./ (1 + exp(-X*b));
    L = y'*log(p) + (1-y')*log(1-p);
    L = -L;
    
    % Gradient  
    dLdb = -(y-p)'*X;
    
end