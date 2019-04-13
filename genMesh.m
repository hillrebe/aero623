function dat = genMesh(maxCamber, locCamber, thickness)

% Inputs:
% 1. maximum camber
% 2. location of maximum camber
% 3. maximum thickness
%
% Outputs:
% 1. structure containing all relevant mesh data

tic
% Settings
plotstuff = 1; %0 to skip plotting mesh, 1 to plot mesh and airfoil
dxLE = 50.; %x-offset of left side of mesh from LE of airfoil
dyhalf = 50.; %y-offset of top/bottom of mesh from airfoil
dxwake = 50.; %x-offset of wake from TE of airfoil
nwake = 100; %Number of nodes in the x-direction of the wake (log spacing)
nfar = 100; %Number of nodes in the y-direction away from the airfoil (log spacing)

% Call genAirfoil to obtain full set of airfoil coordinates
coords = genAirfoil(maxCamber, locCamber, thickness);
xLE = min(coords(:,1)); %x-coordinate of LE nose

% Keep only the odd nodes (80 nodes total)
icoords_keep = 1;
coords_keep = zeros((length(coords)+1)/2,2);
for i = 1:length(coords)
    if mod(i,2) ~= 0
        coords_keep(icoords_keep,:) = coords(i,:);
        icoords_keep = icoords_keep+1;
    end
end
coords = coords_keep;
numpts = length(coords); %Number of airfoil points kept (including double TE, excluding LE)

% Determine edge points - top half

% Apply cosine spacing on curved portion of boundary near LE
numtop = numpts/2 + 1;
angles = linspace(0,pi,numtop)';
xtop = coords(1,1)-((cos(angles)+1.0)/2)*(coords(1,1)+dxLE);
xtop = [xtop(2:end-1);coords(1,1)];
ytop = zeros(length(xtop),1);

for i = 1:length(xtop)
    if xtop(i) >= 0.0
        ytop(i) = dyhalf;
    else
        ytop(i) = sqrt((1-xtop(i)^2/dxLE^2)*dyhalf^2);
    end
end
xtop = flip(xtop); %Right to left
ytop = flip(ytop); %Right to left

% Apply linear spacing on linear portion of boundary
ilin = find(ytop == dyhalf,1,'last'); %Index of last node to be on linear portion
xtoplin = linspace(coords(1,1),xtop(ilin),ilin);
for i = 1:ilin
    xtop(i) = xtoplin(i);
end

% Determine edge points - bottom half

% Apply cosine spacing on curved portion of boundary near LE
numbot = numpts/2 + 1;
angles = linspace(0,pi,numbot)';
xbot = coords(1,1)-((cos(angles)+1.0)/2)*(coords(1,1)+dxLE);
xbot = [xbot(2:end-1);coords(1,1)];
ybot = zeros(length(xbot),1);

for i = 1:length(xbot)
    if xbot(i) >= 0.0
        ybot(i) = -dyhalf;
    else
        ybot(i) = -sqrt((1-xbot(i)^2/dxLE^2)*dyhalf^2);
    end
end

% Apply linear spacing on linear portion of boundary
ilin = find(ybot == -dyhalf,1,'first'); %Index of first node to be on linear portion
xbotlin = linspace(xbot(ilin),coords(1,1),length(xbot)-ilin+1);
for i = ilin:length(xbot)
    xbot(i) = xbotlin(i-ilin+1);
end

%% Wake region

xwake = logspace(log10(coords(1,1)),log10(dxwake),nwake)';
xwake = xwake(2:end);
ywake = coords(1,2)*ones(size(xwake)); %Align wake y-coordinates with TE of airfoil
ywake_far = dyhalf*ones(size(xwake));

%% Assemble full list of coordinates
coords_airfoil = [[flip(xwake),flip(ywake)];coords(1:numpts/2,:);[xLE,0];coords(numpts/2+1:end,:);[xwake,ywake]];
coords_far = [[flip(xwake),flip(ywake_far)];[xtop,ytop];[-dxLE,0];[xbot,ybot];[xwake,-ywake_far]];

% Compute x and y distances of each segment from airfoil/wake centerline to farfield
distx = zeros(length(coords_airfoil),1);
disty = zeros(length(coords_airfoil),1);
for i = 1:length(distx)
    distx(i) = coords_far(i,1)-coords_airfoil(i,1);
end
for i = 1:length(disty)
    disty(i) = coords_far(i,2)-coords_airfoil(i,2);
end

% Apply log spacing from airfoil towards farfield and store node coords
logfrac = logspace(0,2.0,nfar)-1.0;
logfrac = logfrac/logfrac(end);

nodesx = zeros(length(distx),nfar);
for i = 1:length(distx)
    nodesx(i,:) = logfrac*distx(i)+coords_airfoil(i,1);
end
nodesy = zeros(length(disty),nfar);
for i = 1:length(disty)
    nodesy(i,:) = logfrac*disty(i)+coords_airfoil(i,2);
end

% Create V matrix (includes doubled nodes in wake)
V = zeros(size(nodesx,1)*size(nodesx,2),2);
iv = 1;
for i = 1:nfar
    for j = 1:length(distx)
        V(iv,:) = [nodesx(j,i),nodesy(j,i)];
        iv = iv+1;
    end
end

% Delete duplicate nodes from TE through wake
idup = length(distx)-nwake+1:length(distx); %Indexes to remove
V(idup,:) = [];

%% Build E2N matrix

if plotstuff
    figure()
    hold on;
    axis equal;
    axis([-dxLE,dxwake,-dyhalf,dyhalf]);
end

E2N = zeros((length(distx)-1)*(nfar-1),4);
nelem = length(E2N);
%Inner loop
for ie = 1:length(distx)-1
    E2N(ie,:) = [ie,ie+length(distx)-length(idup),ie+length(distx)-length(idup)+1,ie+1];
    if ismember(E2N(ie,1),idup)
        E2N(ie,1) = length(distx)-ie+1;
    end
    if ismember(E2N(ie,4),idup)
        E2N(ie,4) = length(distx)-ie;
    end
    
    if plotstuff
        plotnodes = V(E2N(ie,:),:);
        plotnodes = [plotnodes;plotnodes(1,:)]; %#ok<AGROW>
        plot(plotnodes(:,1),plotnodes(:,2),'k-');
    end
end

%Remaining loops
E2Nresume = [length(distx)-length(idup)+1,2*length(distx)-length(idup)+1,2*length(distx)-length(idup)+2,length(distx)-length(idup)+2];
E2N(ie+1,:) = E2Nresume;

if plotstuff
    plotnodes = V(E2N(ie+1,:),:);
    plotnodes = [plotnodes;plotnodes(1,:)];
    plot(plotnodes(:,1),plotnodes(:,2),'k-');
end

for i = ie+2:(length(distx)-1)*(nfar-1)
    if mod(i-1,length(distx)-1) == 0
        E2N(i,:) = E2N(i-1,:)+2;
    else
        E2N(i,:) = E2N(i-1,:)+1;
    end
    
    if plotstuff
        plotnodes = V(E2N(i,:),:);
        plotnodes = [plotnodes;plotnodes(1,:)]; %#ok<AGROW>
        plot(plotnodes(:,1),plotnodes(:,2),'k-');
    end
end

% Plot airfoil
if plotstuff
    plot(coords(:,1),coords(:,2),'r-','linewidth',1.5)
end

%% Compute Area matrix
Area = zeros(nelem,1);
for ie = 1:nelem
    nodes = V(E2N(ie,:),:);
    Area(ie) = polyarea(nodes(:,1),nodes(:,2));
end

%% Compute I2E/B2E matrixes

% Inner loop
iIE = 1;
iBE = 1;
I2E = zeros((length(distx)-2)*(nfar-1)+(length(distx)-1)*(nfar-1)-nfar-80,4);
B2E = zeros(length(distx)-1+(nfar-1)*2+80,3);

for ie = 1:length(distx)-1
    if ie < length(idup) %Wake interior edges
        I2E(iIE,:) = [ie,4,length(distx)-ie,4];
        iIE = iIE+1;
    end
    if (length(idup) <= ie) && (ie <= length(distx)-length(idup)) %Airfoil boundary edges
        B2E(iBE,:) = [ie,4,1];
        iBE = iBE+1;
    end
    if ie == 1 %First element BE
        B2E(iBE,:) = [ie,1,2];
        iBE = iBE+1;
    end
    if ie == length(distx)-1 %Last element BE
        B2E(iBE,:) = [ie,3,2];
        iBE = iBE+1;
    end
    if ie < length(distx)-1 %Interior edges
        I2E(iIE,:) = [ie,3,ie+1,1];
        iIE = iIE+1;
    end
end

%All other loops
for ie = length(distx):nelem
    
    I2E(iIE,:) = [ie-length(distx)+1,2,ie,4]; %IE (uses previous row)
    iIE = iIE+1;
  
    if mod(ie-1,length(distx)-1) == 0 %First element BE
        B2E(iBE,:) = [ie,1,2];
        iBE = iBE+1;
    end
    
    if mod(ie,length(distx)-1) == 0 %Last element BE
        B2E(iBE,:) = [ie,3,2];
        iBE = iBE+1;
    else
        I2E(iIE,:) = [ie,3,ie+1,1]; %Interior edges
        iIE = iIE+1;
    end
end

%Outer boundary only
for ie = nelem-length(distx)+2:nelem
    B2E(iBE,:) = [ie,2,2];
    iBE = iBE+1;
end

%% Compute In and Bn matrixes
In = zeros(length(I2E),2);
Bn = zeros(length(B2E),2);

for i = 1:length(In)
    elem = I2E(i,1);
    edge = I2E(i,2);
    
    switch edge
        case 1
            inodes = [1,2];
        case 2
            inodes = [2,3];
        case 3
            inodes = [3,4];
        case 4
            inodes = [4,1];
    end
    
    node1 = V(E2N(elem,inodes(1)),:);
    node2 = V(E2N(elem,inodes(2)),:);
    dl = node2-node1;
    n = cross([dl,0],[0,0,1]);
    In(i,:) = n(1:2)'/norm(n); %Normalize       
end

for i = 1:length(Bn)
    elem = B2E(i,1);
    edge = B2E(i,2);
    
    switch edge
        case 1
            inodes = [1,2];
        case 2
            inodes = [2,3];
        case 3
            inodes = [3,4];
        case 4
            inodes = [4,1];
    end
    
    node1 = V(E2N(elem,inodes(1)),:);
    node2 = V(E2N(elem,inodes(2)),:);
    dl = node2-node1;
    n = cross([dl,0],[0,0,1]);
    Bn(i,:) = n(1:2)'/norm(n); %Normalize       
end

%% Store relevant data in struct
dat.V = V;
dat.E2N = E2N;
dat.Area = Area;
dat.I2E = I2E;
dat.B2E = B2E;
dat.In = In;
dat.Bn = Bn;

%% Watertight check
Area_theo = pi*dxLE*dyhalf/2+dxwake*dyhalf*2; %Theoretical (curved) domain area
Area_act = sum(Area); %Actual element areas
Area_pct = Area_act/Area_theo;
if abs(1-Area_pct)>0.0004
    error('ERROR: Mesh may not be watertight - domain area is %.2f percent filled\n',Area_pct*100);
else
    fprintf('Domain area is %.2f percent filled\n',Area_pct*100);
end

fprintf('Mesh generation and plotting completed in %.2f seconds\n',toc);