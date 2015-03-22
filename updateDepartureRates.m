function pathDepartureRates  = updateDepartureRates(N, pathDepartureRates)

% variable 'N' is the number of vehicles on each link at current time step
% i.e. N(1) = Nup(1,t) - Ndown(1,t) 

% random change to pathDepartureRates each time-step
pathDepartureRates = pathDepartureRates + (2*rand(1,1)-1) * 0.25 * pathDepartureRates;


end