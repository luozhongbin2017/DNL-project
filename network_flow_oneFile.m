% Final Year Project - Large Scale Implementation of Dynamic Traffic
% Assignment Algorithms

% Undergraduate Student: Gabriel Eve
% Supervisor: Dr Ke Han

% gabrielfeve@gmail.com
% May 2014

% ======================================================================= %

% Notes and assumptions:
% - Free flow time and equivalent for kinematic wave velocity are multiples
% of dt (if not in input data, then uses modified data in simulation)
% - Sources and sinks represented at OD nodes by means of virtual link
% - Departure rates for all sources
% - Use S.I units i.e. metres, seconds and vehicles

% TODO: are all units correct? do all calculations in S.I.?

clear


%% DATA INPUT

% load network data
uiopen([cd, '/Processed Data/*_FA_pp.mat'])
fprintf(['[1] Loading network data from ', filename, '_FA_pp.mat \n'])

%% Inputs (user changeable)
% TODO: consistency with units

dt = min(Tff);                                % [s] set the time-step
nt = input('Number of time steps: ');          % num of time steps
typDepartureRate = 80;           % [veh/s] source departure rate
typDepartureDur = 30;           % typical number of time steps of departure flow


%% INITIALIZING VARIABLES
fprintf('[2] Initialising variables \n')
if sum(Tff < dt) > 0
    warning('dt > the minimum link free-flow time\ndt = %1.5f, fft_min = %1.5f', dt, min(Tff))
end
% round free flow times to multiples of dt and update L accordingly
Tff_ = round(Tff / dt) * dt;            % [s] modified free flow time
Tff_(Tff_==0) = dt;
L_ = Tff_ .* L ./ Tff;                  % [m] modified length of link
Tkw = 3 * Tff_;                         % [s] (link length) / (kinematic wave velocity)
rho_jL = 4 * C .* Tff_;
rho_j = rho_jL ./ L_;                   % [veh/m] jam density
tnk = round(Tff_ / dt);                 % num of time-steps for free flow to reach head node
tnw = round(Tkw / dt);                  % num of time-steps for shockwave to propagate to tail node
t = 0:dt:(nt-1)*dt;                     % [s] time vector

Nup = zeros(sinkIdx(end), nt);          % [veh] cumulative sum of upstream vehicles
Ndn = zeros(sourceIdx(end), nt);        % [veh] cumulative sum of downstream vehicles
Qin = zeros(sinkIdx(end), nt);          % [veh/s] upstream flow
Qout = zeros(sourceIdx(end), nt);       % [veh/s] downstream flow
% Qinijr = zeros(numTotalPathLinks, nt);  % [veh/s] upstream flow of paths
Nsource = zeros(numSources,1);          % [veh] number of veh queueing at source


%% Departure Rates
fprintf('[3] Setting departure rates \n')

sourceDepartures = zeros(numSources,nt);
startTimes = randi(typDepartureDur,numSources,1);
endTimes = typDepartureDur + randi(typDepartureDur,numSources,1);
depRate = typDepartureRate + (rand(numSources,1)-0.5) .* (typDepartureRate .* 0.25);

for n = 1:numSources
    st = startTimes(n);
    et = endTimes(n);
    sourceDepartures(n,st:et) = depRate(n);
    %     i = sources(n);
    %     sourceDepartures(n,:) = sum(pathDepartures(sourceNode == i,:), 1);  % sum path departure rates for each source
end

Qin(sourceIdx,:) = sourceDepartures;  % set total departure rate from each source

clearvars st et

%% COMPUTATION
% show calculation progress
% h = waitbar(0,'Simulating ...');
fprintf('[4] Simulating traffic flows')
dispstat('', 'init')
dispstat(sprintf('\t0%% Complete'))
CC = [C; inf(numSources,1)];             % flow capacity
SS = [C; inf(numSources+numSinks,1)];    % supply (how much can enter)
S = C;
D = zeros(numLinks,1);
Dsource = zeros(numSources,1);
DD = [D; Dsource];

eps = 1e-9;                         % machine precision tolerance


tic
for tn = 1:nt  % loop over all time steps
    
    tkappa = tn - tnk;                  % t-L/k index
    tomega = tn - tnw;                  % t-L/w index
    
    % demand and suppply on each link and node
    
    Qsource = Qin(sourceIdx,tn);
    
    % Demand capacity produced by link
    ak0 = (tkappa >= 1);
    ak = find(ak0);
    aktk = sub2ind([sinkIdx(end),nt], ak, tkappa(ak0));
    
    KstateD = (Nup(aktk) - Ndn(ak,tn) > eps);  % PDE solution
    
    D(:) = 0;
    D(ak(KstateD)) = C(ak(KstateD));
    D(ak(~KstateD)) = Qin(aktk(~KstateD));
    
    % Supply capacity of link
    aw0 = (tomega >= 1);
    aw = find(aw0);
    awtw = sub2ind([sourceIdx(end),nt], aw, tomega(aw0));
    
    KstateS = (Nup(aw,tn) - Ndn(awtw) - rho_jL(aw) < -eps);  % PDE solution
    
    S(:) = C;
    S(aw(~KstateS)) = Qout(awtw(~KstateS));
    
    
    % determine demand from departure rates at each source
    queue = Nsource > eps;
    Dsource(queue) = inf;
    Dsource(~queue) = Qsource(~queue);
    DD(linkIdx) = D;
    DD(sourceIdx) = Dsource;
    
    % flow in and out of junction
    SS(linkIdx) = S;                         % supply (how much can enter)
    
    for i = 1:numNodes
        nLin = numLinksIn(i);
        Lin = linksIn{i};
        Lout = linksOut{i};
        eta = etas{i};
        
        alpha = alphas{i};
        gamma = cell(nLin,nLout);
        for ik = 1:nLin
            Lik = Lin(ik);
            qin = Qin(Lik,1:tn);
            tau = find(qin(Nup(Lik,1:tn) <= Ndn(Lik,tn) ), 1, 'last');  % logical indexing may be slowing calculation time since not the same size array
            if isempty(tau)
                tau = 1;  % set equal to 1 when there are no values of tau (i.e. before flow has reached end of link)
            end
            qi = Qin(Lik,tau);
            for jk = 1:nLout
                qijr = Qinijr(pathLinksIn{i}{ik,jk}, tau);
                gamma{ik,jk} = max(qijr ./ qi, 0);  % this is because when qi=0 gamma should be 0/is this right? this will mean that the error will always be 'loss' of flow
                alpha(ik,jk) = sum(gamma{ik,jk});
            end
        end

        
        Seff = min([ CC(Lin), SS(Lout,ones(nLin,1))' ./ alpha ],[], 2);
        qout = min(DD(Lin), eta .* Seff);
        Qin(Lout,tn) = alpha' * qout;           % matrix multiplication
        Qout(Lin,tn) = qout;
    end
    
    
    % update number of veh have passed entrance and exit of each link
    Nup(:,tn+1) = Nup(:,tn) + Qin(:,tn) * dt;
    Ndn(:,tn+1) = Ndn(:,tn) + Qout(:,tn) * dt;
    
    % update length of queue at source
    Nsource = max(0, Nsource + dt*(Qin(sourceIdx,tn) - Qout(sourceIdx,tn)));
    
    % update loop progress
    %     waitbar(tn/nt, h)
    dispstat(sprintf('\t%6.0f%% Complete', tn/nt*100))
    
    
end
timeElapsed = toc;
fprintf('\b\tTime elapsed: %6.2fs\n', timeElapsed)

% remove values at Nt+1 for data uniformity
Nup(:,nt+1) = [];
Ndn(:,nt+1) = [];

% clear temp variables
clearvars i n r

% close(h)
% toc

%% OUTPUTS

% average density on each link
rho = (Nup(linkIdx,:) - Ndn(linkIdx,:)) ./ L_(:,ones(1,nt));

%%

fprintf('End \n')

%% debugging
% s_ = pathSinkLinkIdx';
% % Ntotalpathssource = sum(sum(pathDepartures, 2)) * dt
% Ntotalpathssource = sum(sum(Qinijr([r_{:}],:), 2)) * dt
% Ntotalpathsink = sum(sum(Qinijr([s_{:}],:), 2)) * dt
% Ntotalpathsink/Ntotalpathssource
%
% Ntotalsource = sum(sum(Qin(sourceIdx,:), 2)) * dt
% Ntotalsink = sum(sum(Qin(sinkIdx,:), 2)) * dt
% Ntotalsink/Ntotalsource






