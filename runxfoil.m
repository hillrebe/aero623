load('sample_matrix.mat');
nSamples = size(sample_matrix,1);

fid = fopen('/Users/beckyhill/Documents/Aero623/Project4/inputfile.dat', 'wt');

for i = 1:nSamples
    maxCamber = sample_matrix(i,1);
    locCamber = sample_matrix(i,2);
    thickness = sample_matrix(i,3);
    coords = genAirfoil(maxCamber, locCamber, thickness, i);
    filename = sprintf('/Users/beckyhill/Documents/Aero623/Project4/airfoilcoord%.0f.txt \n',i);

    % loading the file into xfoil
    fprintf(fid,'load \n');
    fprintf(fid,filename);

    % moving from xfoil top level to operation level
    fprintf(fid,'OPER \n');

    % for first sample, change from inviscid to viscous, set mach number,
    % set maximum number of iterations
    if i == 1
        fprintf(fid,'ITER 200 \n');  
        fprintf(fid,'M 0.2 \n');
    else    % for all other samples, initialize boundary layers and set reynolds number
        fprintf(fid,'init \n');
    end    

    % record cp and cf data
    fprintf(fid,'a 0 \n');
    fprintf(fid,'cpwr \n');
    fileCp = sprintf('/Users/beckyhill/Documents/Aero623/Project4/cp_%g.dat \n',i);
    fprintf(fid,fileCp);

    % begin polar accumulation
    fprintf(fid,'PACC \n');

    % create file where polar data should be saved
    fileSave = sprintf('/Users/beckyhill/Documents/Aero623/Project4/polar_%g.dat \n',i);
    fprintf(fid,fileSave);
    fprintf(fid, '\n');

    % set angle of attack and generate polar data
    fprintf(fid,'A 0 \n');

    % end polar accumulation
    fprintf(fid, 'PACC \n');

    % delete old polars and return to top level of xfoil
    fprintf(fid, 'pdel 0 \n');
    fprintf(fid, '\n');
end

% close xfoil and close the file
fprintf(fid,'QUIT \n');
fclose(fid);

% run the file through xfoil
cmd = sprintf('cd /Users/beckyhill/Documents/Aero590 && ./xfoil < /Users/beckyhill/Documents/Aero623/Project4/inputfile.dat > /Users/beckyhill/Documents/Aero623/Project4/out.dat');
system(cmd);

% initialize results matrices
polar_results = zeros(nSamples, 3);
cp_results = zeros(80,nSamples);
% cf_results = zeros(159,numSamples);

for j = 1:nSamples

    % open the polar data file that was created for each sample
    fileSaved_pol = sprintf('/Users/beckyhill/Documents/Aero623/Project4/polar_%g.dat',j);
%     fileSaved_cf = sprintf('cf_%g.dat',j);
    fileSaved_cp = sprintf('/Users/beckyhill/Documents/Aero623/Project4/cp_%g.dat',j);

    fid_pol = fopen(fileSaved_pol);
%     fid_cf = fopen(fileSaved_cf);
    fid_cp = fopen(fileSaved_cp);

    % polar data file has extra information not needed in this project, so
    % skip over those lines
    for k = 1:12
        fgetl(fid_pol);
    end

%     fgetl(fid_cf);
    fgetl(fid_cp);

    % extract the coefficients saved in the file
    polar = fscanf(fid_pol, '%g', [7 inf]);
    polar = polar';
    fclose(fid_pol);

%     dump = fscanf(fid_cf, '%g', [8 inf]);
%     cf = dump(7,1:159);
%     cf_results(:,j) = cf';
%     fclose(fid_cf);

    cpdump = fscanf(fid_cp, '%g', [2 inf]);
    cp = cpdump(2,:);
    cp_results(:,j) = cp';
    fclose(fid_cp);

    if j == 1
        x_coord = cpdump(1,:);
    end

    if isempty(polar) == 1 % if the polar is empty the configuration did not converge
        polar_results(j,:) = 0;
        cp_results(:,j) = 0;
%         cf_results(:,j) = 0;
    else % otherwise, add the correct coordinates to the results matrix
        polar_results(j,1) = polar(1,2); %Cl
        polar_results(j,2) = polar(1,3); %Cd
        polar_results(j,3) = polar(1,5); %Cm
%         polar_results(j,4) = min(cf(1:80)); %min Cf
    end
end

% initalize matrix of configurations that did not converge
% removed = [];
% 
% for k = numSamples:-1:1
%     % if the results matrix only has 0s in a row, the configuration did not
%     % converge. delete the row corresponding to this configuration in both
%     % the sample matrix and the results matrix
%     if polar_results(k,1) == 0
%         polar_results(k,:) = [];
%         cf_results(:,k) = [];
%         cp_results(:,k) = [];
%         sample_matrix(k,:) = [];
% 
%         % record which configurations were removed
%         removed = [k, removed];
%     end
% end

% cf_results = cf_results(1:75,:);
% 
% numBasis = 8;
% % [U_cf, S_cf, V_cf] = svd(cf_results);
% % basis_cf = U_cf(:,1:numBasis);
% [U_cp, S_cp, V_cp] = svd(cp_results);
% basis_cp = U_cp(:,1:numBasis);
% 
% numConverged = size(cp_results,2);
% 
% % cf_coeff = zeros(numConverged, numBasis);
% cp_coeff = zeros(numConverged, numBasis);
% 
% for i = 1:numConverged
%     for j = 1:numBasis
% %         cf_coeff(i,j) = dot(cf_results(:,i), basis_cf(:,j));
%         cp_coeff(i,j) = dot(cp_results(:,i), basis_cp(:,j));
%     end  
% end
