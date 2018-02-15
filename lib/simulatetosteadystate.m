function [T,Y,tmax] = simulatetosteadystate(ribo,model,tmax,Y0,hPR,cPR,cXlist,options);

%% simulatetosteadystate
%
%  Function simulates the given circuit in the specifie host. If the
%  circuit has not reached steady state at the end of the simulation then
%  tmax is increased and the circuit is resimulated.
%
%  Steady state is defined as when the largest change in the derivative
%  vector is less than a value specified by options.threshold.
%
%  [T,Y,tmax] = simulatetosteadystate(ribo,model,tmax,Y0,hPR,cPR,cXlist,options);
%
%   Options structure
%   - threshold - maximal value of the derivative dY/dt which is terlerated at steady state (i.e. when is dY/dt = 0) (default = 1)
%   - dt - value to increase tmax if time needs increasing (default = tmax)
%   - timermax - time for which simulatetosteadystate will run before returning failure warning and terminating (default = 10 minutes)
%   - FLAG - print status (default = 0);
%

%% check options struct and create defaults
try options.threshold;
    threshold = options.threshold;
catch
    threshold = 1;
end

try options.dt;
    dt = options.dt;
catch
    dt = tmax;
end

try options.timermax;
    timermax = options.timermax;
catch
    timermax = 600; % 10 minutes
end

try options.FLAG;
    FLAG = options.FLAG;
catch
    FLAG = 0;
end

%% setup
sstest = 0;
timer = tic;

%% specify host model
% convert string to function handles
% string is passed as a 5 character string 'nname' while file names are in the form
% nnameODE so need to add the ODE to the function name during handle creation.
hostODE = str2func([ribo,'ODE']);

%% while loop to find tmax
while sstest == 0;

    % ensure results are not negative
    odeoptions = odeset('NonNegative',[1:length(Y0)]);
   
    % simulate model
    [T,Y] = ode15s( @(T,Y) hostODE(T, Y, model, hPR, cPR, cXlist), [0,tmax], Y0, odeoptions);
    
    % retreive last row of Y vector
    ssY = Y(end,:);
    
    % set population to zero
    ssY(end) = 0;
    
    % evaluate derivative at the end of the simulation
    dY = hostODE([], ssY, model, hPR, cPR, cXlist);
    
    % evaluate if the derivative is approximately zero
    if max(abs(dY)) < threshold;
        sstest = 1;
    else
        sstest = 0;
    end
    
    % if not at steady state then increase tmax
    if sstest == 0;
        if FLAG == 1;
            warning('Not at steady state. Increasing tmax by dt');
        end
        tmax = tmax + dt;
    end
    
    %% failed to find steadystate
    if toc(timer) > timermax;
        sstest = 1;
        if FLAG == 1;
            warning('Failed to find steady state. Returning long term behaviour');
        end
    end
    
end