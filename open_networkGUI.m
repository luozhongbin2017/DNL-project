% openNetworkGUI

if exist('X', 'var') == 1
    
    dispstat('', 'init')
    dispstat('Opening GUI...')
    s = networkGUI(name,filename,X,Y,headNode,tailNode,adjacencyList,rho,rho_j,nt,Nup,Ndn,Qin,Qout,C,numLinks,numNodes,sources,sinks,sourceIdx,sinkIdx,node2source,node2sink,t);
    dispstat('GUI Loaded')
else
    fprintf('Missing node coordinates \n')
end