function sigm = sigmoid(alpha)
    %SIGMOID Custom smooth transition function for aircraft dynamics.
    %
    %   sigm = SIGMOID(alpha) computes a smoothed sigmoid-like value
    %   based on the angle of attack alpha, using parameters from the
    %   aerodynamic model. 
    %
    %   Inputs:
    %       alpha - Angle of attack 
    %
    %   Outputs:
    %       sigm - Smoothed sigmoid-like value

    % Retrieve parameters from the aircraft model
    params = getAircraftParams(); 
    M = params.M;           % Steepness of the transition curve
    alpha0 = params.alpha0; % Central transition angle

    % Compute exponentials for transitions near ±alpha0
    neg_exp = exp(-M * (alpha - alpha0));
    pos_exp = exp(M * (alpha + alpha0));
    
    % Compute the customized sigmoid-like function
    sigm = (1 + neg_exp + pos_exp) / ((1 + neg_exp) * (1 + pos_exp));
end
