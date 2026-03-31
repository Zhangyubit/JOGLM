function [S, I, R, E ] = compartmental_epidemics( A, p, T, seeds, params )
% [S, I, R, E ] = compartmental_epidemics( A, p, T, seeds, params )
%  Simulates an epidemic spread according to the SIR model.
%  A is the graph adjacency matrix,
%  p is the contagion probability,
%  T is the length of the infectious period, and
%  seeds are the vertices infected at time zero (patient zero)

n = size(A, 1);

if not(isfield(params, 'population')), params.population = 1;   end
if not(isfield(params, 'maxTime')),    params.maxTime = round(1000);   end
if not(isfield(params, 'model')),      params.model = 'SIR';   end
if not(isfield(params, 'immunity')),   params.immunity = 1;   end
if not(isfield(params, 'latency')),    params.latency = 1;   end
if not(isfield(params, 'complete')),   params.complete = 0; end

I = zeros(n, params.maxTime+1);
R = zeros(n, params.maxTime+1);
E = zeros(n, params.maxTime+1);

fixedpt = 0;

% scale up
M  = kron(eye(n),ones(1,params.population))/params.population;

if params.complete
    As = kron(A, ones(params.population));
else
    param.connected = 1;
    Gpop = gsp_erdos_renyi(params.population,0.5,param);
    As   = kron(A, Gpop.W);
end

ns = n*params.population;
Is = kron(I, ones(params.population,1)); Is(seeds*params.population,1)=1;
Ss = ~Is;
Rs = kron(R, ones(params.population,1));
Es = kron(E, ones(params.population,1));

% time of infection and recovery
Ts_infect  = nan(ns,1); Ts_infect(logical(Is)) = 0;
Ts_recover = nan(ns,1);
Ts_exposed = nan(ns,1);

Te = params.latency;
Tr = params.immunity;

% run the model
for t = 1:params.maxTime,
    
    switch params.model,
        case 'SI'
            Inew = (As*Is(:,t)*p > rand(ns,1)) .* (Is(:,t)==0) ;
            Rnew = ((t-Ts_infect)>T) .* (Is(:,t)==1) .* (Rs(:,t)==0);
            Enew = 0;
            Snew = 0;
            Ss(:,t+1) = Ss(:,t) - Inew + Rnew;
            Is(:,t+1) = Is(:,t) + Inew - Rnew;
            
        case 'SIR',
            Inew = (As*Is(:,t)*p > rand(ns,1)) .* (Is(:,t)==0) .* (Rs(:,t)==0);
            Rnew = ((t-Ts_infect)>T) .* (Is(:,t)==1) .* (Rs(:,t)==0);
            Enew = 0;
            Snew = 0;
            
            Ss(:,t+1) = Ss(:,t) - Inew;
            Is(:,t+1) = Is(:,t) + Inew - Rnew;
            Rs(:,t+1) = Rs(:,t) + Rnew;
            
        case 'SIRS',
            Inew = (As*Is(:,t)*p > rand(ns,1)) .* (Is(:,t)==0) .* (Rs(:,t)==0);
            Rnew = ((t-Ts_infect)>T) .* (Is(:,t)==1) .* (Rs(:,t)==0);
            Snew = ((t-Ts_recover)>Tr).*(Rs(:,t)==1);
            Enew = 0;
            
            Ss(:,t+1) = Ss(:,t) - Inew + Snew;
            Is(:,t+1) = Is(:,t) + Inew - Rnew;
            Rs(:,t+1) = Rs(:,t) + Rnew - Snew;
            
        case 'SEIRS',
            Enew = (As*Is(:,t)*p > rand(ns,1)) .* (Es(:,t)==0) .* (Rs(:,t)==0).* (Is(:,t)==0);
            Inew = ((t-Ts_exposed)>Te).* (Es(:,t)==1) .* (Is(:,t)==0) .* (Rs(:,t)==0);
            Rnew = ((t-Ts_infect)>T)  .* (Is(:,t)==1) .* (Rs(:,t)==0);
            Snew = ((t-Ts_recover)>Tr).* (Rs(:,t)==1);
            
            Ss(:,t+1) = Ss(:,t) - Enew + Snew;
            Es(:,t+1) = Es(:,t) + Enew - Inew;
            Is(:,t+1) = Is(:,t) + Inew - Rnew;
            Rs(:,t+1) = Rs(:,t) + Rnew - Snew;
    end
    
    Ts_infect(find(Inew))  = t+1;
    Ts_recover(find(Rnew)) = t+1;
    Ts_exposed(find(Enew)) = t+1;
    
    if sum(Ss(:,t+1) + Is(:,t+1) + Rs(:,t+1) + Es(:,t+1)) ~= ns,
        error('invalid state');
    end
    
    if ~any(Inew) && ~any(Snew) && ~any(Enew) && ~any(Rnew)
        fixedpt = fixedpt+1;
    else
        fixedpt = 0;
    end
    
    if fixedpt==(Te+50)
        break
    end
    
%     if max(isnan(Ts_infect)) == 0
%         break;
%     end
end

Ss = Ss(:,1:t);
Es = Es(:,1:t);
Is = Is(:,1:t);
Rs = Rs(:,1:t);

% scale down
S  = M*Ss;
I  = M*Is;
R  = M*Rs;
E  = M*Es;



end

