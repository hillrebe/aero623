% nSamples = 50;
% samples = genSamples(3,nSamples);

load('sample_matrix.mat');
nSamples = size(sample_matrix,1);

% run different sets of data points on different computers
Becky = 1;
Matt = 0;

if Becky == 1
    sample1 = 1;
    lastsample = nSamples/2;
elseif Matt == 1
    sample1 = nSamples/2+1;
    lastsample = nSamples;
else
    fprintf('Not a valid input');
end

cl_fvm = zeros(nSamples/2,1);
cp_fvm = zeros(80,nSamples/2);

for i = sample1:lastsample
    tic;
    fprintf('\n \n \n \n \n CURRENT AIRFOIL: %g \n \n',i);
    clear datafile
    maxCamber = sample_matrix(i,1);
    locCamber = sample_matrix(i,2);
    thickness = sample_matrix(i,3);
    datafile = genMesh(maxCamber, locCamber, thickness, i);
    [umat, ~] = FVM(datafile, 100000, 1, 1);
    [cl_fvm(i), cp_fvm(:,i)] = calcOutputs(umat, datafile, 0);
    savecl = sprintf('cl_all_%.0f.mat',Becky);
    savecp = sprintf('cp_all_%.0f.mat',Becky);
    save(savecl,'cl_fvm')
    save(savecp,'cp_fvm')
    % files end in 0 if Matt's data, 1 if Becky's data
    toc
end