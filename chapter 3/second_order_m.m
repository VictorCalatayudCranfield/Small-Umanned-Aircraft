function second_order_m(block)
    setup(block);
end

function setup(block)
    % Número de puertos
    block.NumInputPorts  = 3;
    block.NumOutputPorts = 1;

    % Dimensiones de entrada y salida
    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;

    for i=1:3
        block.InputPort(i).Dimensions        = 1;
        block.InputPort(i).DirectFeedthrough = true;
    end

    block.OutputPort(1).Dimensions       = 1;

    % Número de estados
    block.NumContStates = 2;

    % Tiempo de muestreo continuo
    block.SampleTimes = [0 0];

    % Métodos (callbacks)
    block.SimStateCompliance = 'DefaultSimState';

    block.RegBlockMethod('InitializeConditions', @InitializeConditions);
    block.RegBlockMethod('Derivatives',          @Derivatives);
    block.RegBlockMethod('Outputs',              @Outputs);
end

function InitializeConditions(block)
    block.ContStates.Data = [0; 0];
end

function Derivatives(block)
    % Aquí defines los parámetros del sistema

    x1 = block.ContStates.Data(1);  % dx
    x2 = block.ContStates.Data(2);  % x
    u     = block.InputPort(1).Data;
    zeta  = block.InputPort(2).Data;
    wn    = block.InputPort(3).Data;

    dx1 = -2*zeta*wn*x1 - wn^2*x2 + u;
    dx2 = x1;

    block.Derivatives.Data = [dx1; dx2];
end

function Outputs(block)
    wn = 5.0; % misma wn
    x2 = block.ContStates.Data(2);
    block.OutputPort(1).Data = wn^2 * x2;
end
