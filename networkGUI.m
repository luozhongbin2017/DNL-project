 function out = networkGUI(name,filename,X,Y,headNode,tailNode,adjacencyList,rho,rho_j,nt,Nup,Ndn,Qin,Qout,C,numLinks,numNodes,sources,sinks,sourceIdx,sinkIdx,node2source,node2sink,t)


%% Callback definitions

    function linkSelection_callback(h,~)
        activeLink = find(hlinks == h);
        
        figure('Numbertitle', 'off', ...
            'Name', ['Link ',num2str(activeLink)])
        hsub = subplot(2,1,1);
        plot(t, Nup(activeLink,:), 'r'); hold on
        set(hsub, 'ygrid', 'on')
        plot(t, Ndn(activeLink,:), 'b')
        box off
        xlabel('t (h)'), ylabel('num veh')
        legend('N_{up}', 'N_{down}', 'location', 'best')
        
        subplot(2,1,2)
        plot(t, Qin(activeLink,:), 'r'); hold on
        plot(t, Qout(activeLink,:), 'b')
        %         plot([t(1), t(end)], C(activeLink), 'g')
        box off
        xlabel('t (h)'), ylabel('flow (veh/h)')
        legend('Q_{in}', 'Q_{out}', 'location', 'best')
        
    end

    function nodeSelection_callback(h,~)
        
        activeNode = find(hnodes == h);
        
        if node2source(activeNode)
            figure('Numbertitle', 'off', ...
                'Name', ['Node ',num2str(activeNode),': Source'])
            hsub = subplot(2,1,1);
            plot(t, Nup(sourceIdx(node2source(activeNode)),:), 'r'); hold on
            set(hsub, 'ygrid', 'on')
            plot(t, Ndn(sourceIdx(node2source(activeNode)),:), 'b')
            box off
            xlabel('t (h)'), ylabel('num veh')
            legend('N_{back of queue}', 'N_{entry to network}', 'location', 'best')
            
            subplot(2,1,2)
            plot(t, Qin(sourceIdx(node2source(activeNode)),:), 'r'); hold on
            plot(t, Qout(sourceIdx(node2source(activeNode)),:), 'b')
            box off
            xlabel('t (h)'), ylabel('flow (veh/h)')
            legend('Flow to queue', 'Flow into network', 'location', 'best')
        end
        
        if node2sink(activeNode)
            figure('Numbertitle', 'off', ...
                'Name', ['Node ',num2str(activeNode),': Sink'])
            hsub = subplot(2,1,1);
            plot(t, Nup(sinkIdx(node2sink(activeNode)),:), 'b'); hold on
            set(hsub, 'ygrid', 'on')
            box off
            xlabel('t (h)'), ylabel('num veh')
            legend('N_{exit network}', 'location', 'best')
            
            subplot(2,1,2)
            plot(t, Qin(sinkIdx(node2sink(activeNode)),:), 'b'); hold on
            box off
            xlabel('t (h)'), ylabel('flow (veh/h)')
            legend('Flow out of network', 'location', 'best')
        end
        
    end

    function nodes_callback(h,~)
        
        if strcmp(get(h, 'String'), 'Hide Nodes')
            set(hnodes, 'Visible', 'off')
            set(h, 'String', 'Show Nodes')
        else
            set(hnodes, 'Visible', 'on')
            set(h, 'String', 'Hide Nodes')
        end
        
    end

    function slider_callback(h,~)
        
        set(hplay, 'String', 'Start')
        tval = round(get(h, 'Value'));
        update_ui(tval)
        
    end

    function rewind_callback(~,~)
        
        set(hplay, 'String', 'Start')
        tval = 1;
        update_ui(tval)
        
    end

    function popup_callback(h,~)
        
        ind = get(h, 'Value');
        cmaplabel = popupInfo{2,ind};
        cmap = popupInfo{3,ind};
        cind = popupInfo{4,ind};
        caxis(popupInfo{5,ind})
        colormap(cmap)
        hcol = colorbar;
        ylabel(hcol, cmaplabel)
        if strcmp(get(hplay, 'String'), 'Start')
            update_ui(tval)
        end
        
    end

    function play_callback(h,~)  % *** + GIF generator ***
        
        if strcmp(get(h, 'String'), 'Start')
            set(h, 'String', 'Stop')
            tval = get(hslider, 'Value');
            step = 5;
            
            for tn = tval:step:nt
                if strcmp(get(h, 'String'), 'Start')
                    return
                end
                update_ui(tn)
%                 %%% start save to GIF %%%
%                 drawnow
%                 frame = getframe(hfig);
%                 im = frame2im(frame);
%                 [A,map] = rgb2ind(im,256);
%                 if tn == 1;
%                     imwrite(A,map,[filename, '.gif'],'gif','LoopCount',Inf,'DelayTime',0.2);
%                 else
%                     imwrite(A,map,[filename, '.gif'],'gif','WriteMode','append','DelayTime',0.2);
%                 end
%                 %%% end save to GIF %%%
                pause(0.1)
            end
            set(h, 'String', 'Start')
            
        else  % if the button string is not 'Start'
            set(h, 'String', 'Start')
        end
        
    end

    function update_ui(tn)
        
        set(hslider, 'Value', tn)
        set(hlinks, {'Color'}, num2cell(cmap(cind(:,tn),:), 2))
        set(hnodes, {'MarkerFaceColor'}, num2cell(cmapN(cindN(:,tn),:), 2))
%         sel = sel1(:,tn)|sel2(:,tn);
        %         set(hnodes(sel), 'MarkerEdge', 'none')
        %         set(hnodes(~sel), 'MarkerEdge', 'k', 'MarkerFace', 'none')
%         set(hnodes(~sel), 'MarkerFace', 'none')
        textsliderstr = sprintf('t = %0.2f h', t(tn));
        set(htextslider, 'String', textsliderstr)
        
    end

    function zoom_callback(~,~)
        nscale = get(ax, 'Xlim');
        nscale = nscale(2) - nscale(1);
        zf = scale/nscale;
        set(hlinks, 'LineWidth', zf^0.6*linkSize)
        set(hnodes, 'MarkerSize', zf^0.6 * nodeSize)
        set(hlinks, {'Xdata'}, num2cell(X(adjacencyList) + dx(:,[1 1]) / zf^0.6, 2), ...
            {'Ydata'}, num2cell(Y(adjacencyList) + dy(:,[1 1]) / zf^0.6, 2))
    end


%% Create main GUI figure

networkScale = sqrt(numNodes);
linkSize = 38 / networkScale;  % link width
nodeSize = 125 / networkScale;  % node diameter

scrsz = get(0,'ScreenSize');
hfig = figure('Name', 'Network Simulation', ...
    'OuterPosition',[5 5 scrsz(3)*0.7 scrsz(4)*0.95], ...
    'NumberTitle', 'off');
hpanelmain = uipanel(hfig, 'Position', [0.15 0 0.85 1], 'Background', 'w');
hpanelside = uipanel(hfig, 'Position', [0 0 0.15 1]);
ax = axes('Parent', hpanelmain, 'Position', [0.05 0.05 0.9 0.9]);
title([name, ' Network'], 'FontSize', 22)
hold on

% create virtual nodal coordinates for 2-way flow (left hand drive - for
% right hand drive, reverse the sign for dx and dy.
laneSpacing = (max(X) - min(X)) / networkScale * 0.05;
phi = atan2(Y(headNode) - Y(tailNode), X(headNode) - X(tailNode));
dx = -sin(phi) * laneSpacing;
dy = cos(phi) * laneSpacing;


%% create colourmap for DTA GUI
% uncongested = [linspace(0.7,0,16)', linspace(1,0.7,16)', linspace(0.8,0,16)'];
% congested = [linspace(1,0.8,48)', linspace(1,0,48)', linspace(0,0,48)'];
% cmaprhor  = [uncongested; congested];
cmaprhor = [linspace(0.7,0,16)', linspace(1,1,16)', linspace(0.8,0,16)'; ...
    linspace(0,1,16)', linspace(1,1,16)', linspace(0,0,16)'; ...
    linspace(1,0.8,32)', linspace(1,0,32)', linspace(0,0,32)'];
uncongested = [linspace(0.7,1,16)', linspace(1,1,16)', linspace(0.8,0,16)'];
congested = [linspace(1,0.8,48)', linspace(1,0,48)', linspace(0,0,48)'];
cmaprho  = [uncongested; congested];
% cmaprho = [linspace(0.2,1,32)',linspace(1,1,32)',linspace(0.2,0,32)'; ...
%     linspace(1,0.8,32)',linspace(1,0,32)',linspace(0,0,32)'];
cmapQ = [linspace(0.9,0,32)',linspace(0.9,0,32)',linspace(0.9,1,32)'; ...
    linspace(0,0,32)',linspace(0,0,32)',linspace(1,0,32)'];
cmapN = colormap([0.8 0.8 0.8; 1 0.2 0.2; 0.2 0.8 0.8; 1 0.2 1]);
cmap = [];

maxrho = max(rho(:))/3;
cindrho = min(max(round(size(cmaprho,1) * rho / maxrho), 1), size(cmaprho,1));
cindrhor = min(max(round(size(cmaprhor,1) * rho ./ rho_j(:,ones(1,nt))), 1), size(cmaprhor,1)); % >1 & <100
cindQin = max(round(size(cmapQ,1) * Qin(1:numLinks,:) ./ C(:,ones(1,nt))), 1);  % at least as large as 1
cindQout = max(round(size(cmapQ,1) * Qout(1:numLinks,:) ./ C(:,ones(1,nt))), 1);
% colour coding for nodes
sel1 = false(numNodes,nt); sel2 = sel1;
sel1(sources,:) = Qout(sourceIdx,:) > 0;
sel2(sinks,:) = Qin(sinkIdx,:) > 0;
cindN = ones(size(sel1));
cindN(sel1) = 2; cindN(sel2) = 3; cindN(sel1 & sel2) = 4;
cind = [];

cmaplabelQ = 'Relative Flow (Q/C %)';
cmaplabelrhor = 'Relative Density (\rho/\rho_j %)';
cmaplabelrho = 'Density (\rho veh/mi)';

ratio = [0 100];
tval = 1;

%% Plot Network
% plot links
hlinks = plot(ax, (X(adjacencyList) + dx(:,[1 1]))', ...
    (Y(adjacencyList) + dy(:,[1 1]))');
set(hlinks, 'LineWidth', linkSize)

% plot nodes
hnodes = NaN(numNodes,1);
for n = 1:numNodes
    hnodes(n) = line(X(n), Y(n));
end

% nb = 1;
% ne = min(numNodes,100);
% while nb <= numNodes
%     get(scatter(ax, X(nb:ne), Y(nb:ne)), 'Children')
%     hnodes(nb:ne) = get(scatter(ax, X(nb:ne), Y(nb:ne)), 'Children');
%     nb = ne+1;
%     ne = min(nb+99, numNodes);
% end


hnodes = sort(hnodes);
set(hnodes, 'Marker', 'o', 'MarkerEdge', 'none', 'HitTest', 'on', 'MarkerFaceColor', [0.2 0.2 0.2])

% htext = text(X,Y, cellstr(num2str((1:numNodes)')), 'Horizontal', 'left');
% set(htext, 'HitTest', 'off')

% attach callback functions to nodes and link handles
set(hlinks, 'ButtonDownFcn', @linkSelection_callback)
set(hnodes, 'ButtonDownFcn', @nodeSelection_callback)

set(zoom(hfig), 'ActionPostCallback', @zoom_callback, 'Enable', 'on')
scale = get(ax, 'xlim');
scale = scale(2) - scale(1);

%% UI elements
buttonheight = 0.05;
popupInfo(:,1) = {'Absolute Density'; cmaplabelrho; cmaprho; cindrho; [0 maxrho]};
popupInfo(:,2) = {'Relative Density'; cmaplabelrhor; cmaprhor; cindrhor; ratio};
popupInfo(:,3) = {'Relative flow in'; cmaplabelQ; cmapQ; cindQin; ratio};
popupInfo(:,4) = {'Relative flow out'; cmaplabelQ; cmapQ; cindQout; ratio};

% make a slider
hslider = uicontrol(hpanelmain, 'Style', 'slider', ...
    'Position', [100, 5, 500, 20], ...
    'Min', 1, 'Max', nt, ...
    'Value', tval, ...
    'SliderStep', [0.05 0.05], ...
    'Callback', @slider_callback);

% make text display above slider
htextslider = uicontrol(hpanelmain, 'Style', 'text', ...
    'Background', 'w', ...
    'Position', [5, 7, 90, 20], ...
    'FontSize', 16, ...
    'String', ['t = ', num2str(tval), ' h']);

% make a start/stop button
hplay = uicontrol(hpanelside, 'Style', 'pushbutton', ...
    'Units', 'normal', ...
    'Position', [0.1, 0.85, 0.8, 0.1], ...
    'FontSize', 16, ...
    'String', 'Start', ...
    'Callback', @play_callback);

% make a rewind button
uicontrol(hpanelside, 'Style', 'pushbutton', ...
    'Units', 'normal', ...
    'Position', [0.1, 0.8, 0.8, buttonheight], ...
    'FontSize', 12, ...
    'String', '|<<', ...
    'Callback', @rewind_callback);

% make a rewind button
uicontrol(hpanelside, 'Style', 'pushbutton', ...
    'Units', 'normal', ...
    'Position', [0.1, 0.6, 0.8, buttonheight], ...
    'FontSize', 12, ...
    'String', 'Hide Nodes', ...
    'Callback', @nodes_callback);

% pop-up choice of density, flow
defpop = 2;
hpop = uicontrol(hpanelside, 'Style', 'popupmenu', ...
    'Units', 'normal', ...
    'Position', [0.05, 0.7, 0.9, buttonheight], ...
    'FontSize', 14, ...
    'String', popupInfo(defpop,:), ...
    'Value', defpop, ...
    'Callback', @popup_callback);


%% update UI initially
popup_callback(hpop,[])
zoom_callback([],[])
axis off equal

out = [];
end
