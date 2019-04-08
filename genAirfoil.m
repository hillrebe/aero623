function airfoilCoord = genAirfoil(maxCamber, locCamber, thickness,def)

% Inputs:
% 1. maximum camber
% 2. location of maximum camber
% 3. maximum thickness
% 4. number used to define the airfoil in a text file

% converting the inputs to fractions
maxCamber = maxCamber/100; 
locCamber = locCamber/10;
thickness = thickness/100;

% x coordinates taken from xfoil's NACA 0012 airfoil
x = [1.000000 0.9916770 0.9803640 0.9672594 0.9527094 0.9371906 ...
    0.9211098 0.9047307 0.8881990 0.8715890 0.8549374 0.8382620 0.8215722 ...
    0.8048727 0.7881664 0.7714554 0.7547412 0.7380253 0.7213088 0.7045932 ...
    0.6878796 0.6711694 0.6544637 0.6377640 0.6210716 0.6043879 0.5877143 ...
    0.5710523 0.5544034 0.5377692 0.5211515 0.5045520 0.4879725 0.4714152 ...
    0.4548821 0.4383756 0.4218982 0.4054527 0.3890420 0.3726696 0.3563391 ...
    0.3400549 0.3238215 0.3076446 0.2915302 0.2754859 0.2595202 0.2436435 ...
    0.2278683 0.2122105 0.1966898 0.1813320 0.1661718 0.1512560 0.1366496 ...
    0.1224422 0.1087546 0.9573936e-1 0.8357044e-1 0.7241381e-1 0.6238752e-1 ...
    0.5353091e-1 0.4580157e-1 0.3909828e-1 0.3329271e-1 0.2825461e-1 ...
    0.2386618e-1 0.2002762e-1 0.1665766e-1 0.1369187e-1 0.1108041e-1 ...
    0.8785821e-2 0.6781149e-2 0.5048426e-2 0.3577154e-2 0.2362688e-2 ...
    0.1403607e-2 0.6990502e-3 0.2434807e-3 0.2600131e-4];
   

%thickness distribution
yt = thickness .* (1.4845.*x.^0.5 - 0.63.*x - 1.758.*x.^2 + 1.4215.*x.^3 - 0.5075.*x.^4);

% camber line, before max camber
yc1 = (maxCamber/locCamber^2) .* (2.*locCamber.*x - (x.^2));
dyc1 = 2*maxCamber/locCamber^2 .* (locCamber - x);

% camber line, after max camber
yc2 = (maxCamber / (1-locCamber)^2) .* (1 - 2*locCamber + 2*locCamber.*x - x.^2);
dyc2 = (2*maxCamber / (1-locCamber)^2) .* (locCamber - x);

% initialize coordinate matrix
airfoilCoord = zeros(length(x)*2, 2);

for i = 1:length(x)
    
    % set equations for camber line for before and after loc of max camber
    if x(i) < locCamber
        yc = yc1(i);
        dyc = dyc1(i);
    else
        yc = yc2(i);
        dyc = dyc2(i);
    end
    
    % calculate theta
    theta = atan(dyc);
    
    % upper coordinates
    x_upper = x(i) - yt(i)*sin(theta);
    y_upper = yc + yt(i)*cos(theta);
    
    %lower coordinates
    x_lower = x(i) + yt(i)*sin(theta);
    y_lower = yc - yt(i)*cos(theta);
    
    % add upper coordinates to matrix
    airfoilCoord(i,1) = x_upper;
    airfoilCoord(i,2) = y_upper;
    
    % add lower coordinates to matrix
    airfoilCoord(end+1-i, 1) = x_lower;
    airfoilCoord(end+1-i, 2) = y_lower;
    
end

% airfoilCoord = airfoilCoord./(max(airfoilCoord))
airfoilCoord(length(x),:) = [];

% create a file for the airfoil and populate it
filename = sprintf('airfoilcoord%.0f.txt',def);
fid = fopen(filename, 'w+t');
fprintf(fid,'testfoil \n'); % at the top of the file, name xfoil will use
fprintf(fid,'     %.8f    %.8f \n',airfoilCoord.'); % adding coordinates
fclose(fid);

end
