% Dynamic network loading
% - Script for bundling network data into a mat-file in specified
% format for use with the dynamic network loading (DNL) program

% Fixed-alpha version!!


%% List of inputs


clear
close all

% open data file containing:
% (1) linkData (mx5) [start node, end node, capacity, length, free flow time]
% m is the link number
% (2) [pathList]
% (3*) nodeCoordinates 
% (4*) name (of network)
% (5*) filename
uiopen([cd, '/Raw Data/*_raw.mat'])

sourceEta = 0.3;   % TODO: currently this is a made up quantity
warning(['eta = ', num2str(sourceEta), ' for every source'])

%% Adjacency list/matrix

adjacencyList = linkData(:,1:2);
tailNode = adjacencyList(:,1);
headNode = adjacencyList(:,2);
numLinks = size(linkData, 1);
numNodes = length(unique(adjacencyList(:)));
if numNodes ~= max(adjacencyList(:));
    warning('Nodes either do not start at ''1'' or are not sequential')
end

% obtain matrix from adjacency list
adjacencyMatrix = sparse(tailNode, headNode, true, numNodes, numNodes);


%% Link properties

C = linkData(:,3);                  % [veh/h] link capacity
L = linkData(:,4);                  % [mi] original length of link
Tff = linkData(:,5);                % [h] original free flow time


%% Nodal coordinates (* if available)

if exist('nodeCoordinates', 'var')
X = nodeCoordinates(:,1);
Y = nodeCoordinates(:,2);
end

%% Etas (+ fixed-alphas)

% create sets of ingoing and outgoing links for each node and eta values
% based on ratio of link capacities

% [~,c] = find(pathList');  % list path index for each path item
% c_ = find([diff(c); 1]);  % find index at end of each path
% pathLength = [c_(1); diff(c_)];  % find number of links in each path

% sourceNode = tailNode(pathList(:,1));
% sinkNode = headNode(pathList(sub2ind(size(pathList),1:size(pathList,1),pathLength')));
% sources = unique(sourceNode);
% sinks = unique(sinkNode);
sources = unique(tailNode);
sinks = unique(headNode);
numSources = numel(sources);
numSinks = numel(sinks);
% numPaths = size(pathList,1);


linkIdx = 1:numLinks;
sourceIdx = numLinks + (1:numSources);
sinkIdx = sourceIdx(end) + (1:numSinks);

etas = cell(numNodes,1);
alphas = etas;                 % pre set turning ratios
linksIn = cell(numNodes,1);             % list of links into each node (including sources)
linksOut = cell(numNodes,1);            % list of links out of each node (including sinks)
numLinksIn = zeros(numNodes,1);         % num flows into node
numLinksOut = zeros(numNodes,1);        % num flows out of node
linkInIdx = zeros(sourceIdx(end),1);
linkOutIdx = zeros(sinkIdx(end),1);

B = [headNode; sources];
A = [tailNode; zeros(size(sources)); sinks];
for i = 1:numNodes
    Lin = find(B == i);
    Lout = find(A == i);
    linksIn{i} = Lin;
    linksOut{i} = Lout;
    
    nLin = numel(Lin);
    nLout = numel(Lout);
    numLinksIn(i) = nLin;
    numLinksOut(i) = nLout;
    
    linkInIdx(Lin) = 1:nLin;
    linkOutIdx(Lout) = 1:nLout;
    
    % alphas
    alphas{i} = repmat(1/nLout, nLin, nLout);  % basic alpha generator (all equal for every junction)
    
    if sum(sources == i) == 1            % if node i is a source
        Cin = C(Lin(1:nLin-1));
        etas{i} = [Cin / sum(Cin) * (1-sourceEta); sourceEta];
        
    else
        Cin = C(Lin);
        etas{i} = Cin / sum(Cin);
    end
end

% clear temp variables
clearvars A B Lin Lout nLin nLout Cin

node2source = zeros(numNodes,1);
node2sink = zeros(numNodes,1);
node2source(sources) = 1:numSources;  % sparse?
node2sink(sinks) = 1:numSinks;

%% Paths

%{

% create individual index for path flow on every link
% pre-allocate
sourcePaths = cell(numSources,1);
numPathLinks = numel(find(pathList));
pathSourceLinkIdx = numPathLinks + (1:numPaths);
pathSinkLinkIdx = pathSourceLinkIdx(end) + (1:numPaths);
ipl = 0;  % initialise pathlink index counter


sourceInIdx = linkInIdx(sourceIdx(node2source(sourceNode))); % link-in index of each path source
sinkOutIdx = linkOutIdx(sinkIdx(node2sink(sinkNode)));  % link-out index of each path sink

% process paths
pathLinksIn = cell(numNodes,1);
for i = 1:numNodes
    pathLinksIn{i} = cell(numLinksIn(i), numLinksOut(i));
end
pathLinksOut = pathLinksIn;

for r = 1:numPaths
    pl = pathLength(r);
    path = pathList(r,1:pl);
    kin = sourceInIdx(r);
    kout = linkOutIdx(path(1));
    pathLinksIn{sourceNode(r)}{kin,kout}(end+1) = pathSourceLinkIdx(r);
    ipl = ipl + 1;
    pathLinksOut{sourceNode(r)}{kin,kout}(end+1) = ipl;
    for i = 1:pl-1 
        Lin = path(i);
        Lout = path(i+1);
        kin = linkInIdx(Lin);
        kout = linkOutIdx(Lout);
        pathLinksIn{headNode(Lin)}{kin,kout}(end+1) = ipl;
        ipl = ipl + 1;
        pathLinksOut{headNode(Lin)}{kin,kout}(end+1) = ipl;
    end
    kin = linkInIdx(path(pl));
    kout = sinkOutIdx(r);
    pathLinksIn{sinkNode(r)}{kin,kout}(end+1) = ipl;
    pathLinksOut{sinkNode(r)}{kin,kout}(end+1) = pathSinkLinkIdx(r);
end

numTotalPathLinks = pathSinkLinkIdx(end);

% clear temp variables
clearvars c c_ i ipl kin kout Lin Lout path pl r 

%}

%% Save to file

save([cd, '/Processed Data/', filename, '_FA_pp.mat'])

fprintf(['Variables saved to .../Processed Data/', filename, '_FA_pp.mat', '\n'])
