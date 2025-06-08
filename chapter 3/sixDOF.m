function sixDOF(block)
    setup(block);
end

function setup(block)
    % Número de puertos

    block.NumInputPorts  = 6;
    block.NumOutputPorts = 2;

    % Dimensiones de entrada y salida
    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;

    for i=1: block.NumInputPorts
        block.InputPort(i).Dimensions        = 1;
        block.InputPort(i).DirectFeedthrough = true;
    end
    
    
    block.OutputPort(1).Dimensions  = 3;
    block.OutputPort(2).Dimensions  = 3;
    

    % Número de estados
    block.NumContStates = 13;

    % Tiempo de muestreo continuo
    block.SampleTimes = [0 0];

    % Métodos (callbacks)
    block.SimStateCompliance = 'DefaultSimState';

    block.RegBlockMethod('InitializeConditions', @InitializeConditions);
    block.RegBlockMethod('Derivatives',          @Derivatives);
    block.RegBlockMethod('Outputs',              @Outputs);
end

function InitializeConditions(block)
    block.ContStates.Data = zeros(13,1);
    block.ContStates.Data(7:10)=eul2quat([0,0,0]);
end

function Derivatives(block)

    % external forces
    fx = block.InputPort(1).Data; 
    fy = block.InputPort(2).Data; 
    fz = block.InputPort(3).Data; 
    % external moments
    l = block.InputPort(4).Data;
    m = block.InputPort(5).Data;
    n = block.InputPort(6).Data;
    
    % aircraft parameters
    mass  = 1.56; % mass
    jx = 0.1147; % inertia
    jy = 0.0576; % inertia
    jz = 0.1712; % inertia
    jxz= 0.0015; % inertia


    pn = block.ContStates.Data(1);  % position north
    pe = block.ContStates.Data(2);  % position east
    pd = block.ContStates.Data(3);  % position down

    u = block.ContStates.Data(4);  % u vel 
    v = block.ContStates.Data(5);  % v vel
    w = block.ContStates.Data(6);  % w vel

    e0 = block.ContStates.Data(7);   % quat
    e1 = block.ContStates.Data(8);   % quat
    e2 = block.ContStates.Data(9);   % quat
    e3 = block.ContStates.Data(10);  % quat

    p = block.ContStates.Data(11);  % rot vel
    q = block.ContStates.Data(12);  % rot vel
    r = block.ContStates.Data(13);  % rot vel


    % precompute inertial terms for readability
    gamma  = jx*jz -jxz^2;
    gamma1 = jxz*(jx -jy + jz)/gamma;
    gamma2 = (jz*(jz -jy) + jxz^2)/gamma;
    gamma3 = jz/gamma;
    gamma4 = jxz/gamma;
    gamma5 = (jz -jx)/jy;
    gamma6 = jxz/jy;
    gamma7 = ((jx -jy)*jx + jxz^2)/gamma;
    gamma8 = jx/gamma;
    
    % to ensure norm e =1
    norme=norm([e0,e1,e2,e3]);
    lambda=1000;

    % 12 state model
    diffp = [e1^2+e0^2-e2^2-e3^2, 2*(e1*e2 -e3*e0), 2*(e1*e3+e2*e0);...
             2*(e1*e2+e3*e0), e2^2+e0^2-e1^2-e3^2, 2*(e2*e3-e1*e0);...
             2*(e1*e2-e3*e0), 2*(e2*e3+e1*e0), e3^2+e0^2-e1^2-e2^2]*[u, v, w]';
    
    diffvel = [r*v-q*w; p*w-r*u; q*u-p*v]+1/mass*[fx;fy;fz];
    
    diffquat = [lambda*(1-norme^2), -p, -q, -r; p, lambda*(1-norme^2), r, -q; q, -r, lambda*(1-norme^2), p; r, q, -p, lambda*(1-norme^2)]*0.5*[e0; e1; e2; e3]; 
    
    diffrotvel = [gamma1*p*q-gamma2*q*r; gamma5*p*r-gamma6*(p^2-r^2); gamma7*p*q-gamma1*q*r] +...
                 [gamma3*l+gamma4*n ; m/jy; gamma4*l+gamma8*n];


    block.Derivatives.Data = [diffp; diffvel;diffquat;diffrotvel];
    disp("uDot")
    disp(diffvel(1))
    disp("u")
    disp(u)
    disp("pn")
    disp(pn)
end

function Outputs(block)
    if any(isnan(block.ContStates.Data))
    disp("WARNING: NaN detected in states!");
    end
    block.OutputPort(1).Data = block.ContStates.Data(1:3);
    eul=quat2eul(quaternion(block.ContStates.Data(7:10)'),"XYZ");
    block.OutputPort(2).Data = eul;
    %disp("outputs")
    %disp(block.ContStates.Data(4:6))

    
end
