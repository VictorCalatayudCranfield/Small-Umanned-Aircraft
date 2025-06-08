% forces_moments.m
%   Computes the forces and moments acting on the airframe. 
%
%   Output is
%       F     - forces
%       M     - moments
%       Va    - airspeed
%       alpha - angle of attack
%       beta  - sideslip angle
%       wind  - wind vector in the inertial frame
%

function out = forces_moments(x, delta, wind, params)

%% Input settings    
% relabel the inputs
    pn      = x(1);
    pe      = x(2);
    pd      = x(3);
    u       = x(4);
    v       = x(5);
    w       = x(6);
    phi     = x(7);
    theta   = x(8);
    psi     = x(9);
    p       = x(10);
    q       = x(11);
    r       = x(12);
    delta_e = delta(1);
    delta_a = delta(2);
    delta_r = delta(3);
    delta_t = delta(4);
    w_ns    = wind(1); % steady wind - North
    w_es    = wind(2); % steady wind - East
    w_ds    = wind(3); % steady wind - Down
    u_wg    = wind(4); % gust along body x-axis
    v_wg    = wind(5); % gust along body y-axis    
    w_wg    = wind(6); % gust along body z-axis

    % load aircraft parameters
    struct2vars(params)



    %% Gravity Force component
    Fgrav = [-mass*grav*sin(theta);...
            mass*grav*cos(theta)*sin(phi);...
            mass*grav*cos(theta)*cos(phi)];
    %% Aerodinamic Force Components
    % Aerodynamic coefficients

        %%% calculate alpha
        % 1 Build the quaternion for body→NED from ZYX Euler angles 
        q_b2n = quaternion([psi, theta, phi], 'eulerd', 'ZYX', 'frame');

        % 2 Invert it to get NED→body
        q_n2b = conj(q_b2n);

        % 3 Rotate the steady wind vector (NED) into the body frame
        w_ned = [w_ns, w_es, w_ds];      % row vector
        wind_b_steady = rotatepoint(q_n2b, w_ned);

        % 4 Add gust (already in body axes)
        gust_b = [u_wg, v_wg, w_wg];
        V_wind_b = wind_b_steady + gust_b;  

        % 5 Compute body-relative air velocity
        V_air_vec  = [u, v, w] - V_wind_b;   
        V_air=norm(V_air_vec);

        % 6 Angle of attack
        alpha = atan2(V_air_vec(3), V_air_vec(1));
        beta  = asin(V_air_vec(2)/V_air);

    CD = CD_0 + (CL_0 + CL_alpha*alpha)^2/(pi*e*AR);
    CL = (1-sigmoid(alpha))*(CL_0+CL_alpha*alpha)+sigmoid(alpha)*(2*sign(alpha)*sin(alpha)^2*cos(alpha));
    
    CX         = -CD        *cos(alpha) + CL        *sin(alpha);
    CX_q       = -CD_q      *cos(alpha) + CL_q      *sin(alpha);
    CX_delta_e = -CD_delta_e*cos(alpha) + CL_delta_e*sin(alpha);
    CZ         = -CD        *sin(alpha) - CL        *cos(alpha);
    CZ_q       = -CD_q      *sin(alpha) - CL_q      *cos(alpha);
    CZ_delta_e = -CD_delta_e*sin(alpha) - CL_delta_e*cos(alpha);

    % Aerodinamic Force Components
    Faero = 0.5*rho*V_air^2*S*...
        [CX+CX_q*c/(2*V_air)*q+CX_delta_e*delta_e;...
         CY_0+CY_beta*beta+CY_p*b/(2*V_air)*p+CY_r*b/(2*V_air)*r+CY_delta_a*delta_a+CY_delta_r*delta_r;... 
         CZ+CZ_q*c/(2*V_air)*q+CZ_delta_e*delta_e];

    %% Propeler forces
    Fprop = 0.5*rho*S_prop*Cprop*[(k_motor*delta_t)^2-V_air^2;0;0];
    
    %% Total Force
    Force = Fgrav + Faero + Fprop;
    
    %% Aerodinamic Moments

    Maero = 0.5*rho*V_air^2*S*...
        [b*(Cl_0+Cl_beta*beta+Cl_p*b/(2*V_air)*p+Cl_r*b/(2*V_air)*r+Cl_delta_a*delta_a+Cl_delta_r*delta_r);...
         c*(Cm_0+Cm_alpha*alpha+Cm_q*c/(2*V_air)*q+Cm_delta_e*delta_e);...
         b*(Cn_0+Cn_beta*beta+Cn_p*b/(2*V_air)*p+Cn_r*b/(2*V_air)*r+Cn_delta_a*delta_a+Cn_delta_r*delta_r)];

    %% Propeler Moments

    Mprop=[-k_T_P*(k_Omega*delta_t)^2;0;0];

    %% Total Moments
    Torque = Maero + Mprop;
   
    out = [Force; Torque; V_air; alpha; beta; w_ned(1); w_ned(2); w_ned(3)];
end



