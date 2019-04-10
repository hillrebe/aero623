function dat = genMesh(maxCamber, locCamber, thickness)

% Inputs:
% 1. maximum camber
% 2. location of maximum camber
% 3. maximum thickness
%
% Outputs:
% 1. structure containing all relevant mesh data

% Settings
dxLE = 1; %x-offset of left side of mesh from LE of airfoil
dyhalf = 1; %y-offset of top/bottom of mesh from airfoil

% Call genAirfoil to obtain full set of airfoil coordinates
coords = genAirfoil(maxCamber, locCamber, thickness);

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

% Determine edge points - top half
numtop = icoords_keep/2 + 1;
angles = linspace(0,pi,numtop)';
xtop = 1.0-((cos(angles)+1.0)/2)*(1.0+dxLE);
xtop = [xtop(2:end-1);1.0];
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

% Plot for debugging
figure()
plot(coords(:,1),coords(:,2),'r-','linewidth',1.5)
hold on;
axis equal;
axis([-1,2,-1,1]);
for i = 1:icoords_keep/2
    plot([coords(i,1),xtop(i)],[coords(i,2),ytop(i)],'k-')
end

plot(xtop,ytop,'k.')

% Store relevant data in struct
dat.V = coords;
