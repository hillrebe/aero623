function sample_matrix = genSamples(numParams, samples)

%     global sample_matrix
    
% parameters always in the same order:
%   1. max camber (0-6)
%   2. location of max camber (2-6)
%   3. max thickness (6-12)
%   4. Reynolds number (5-7)
%   5. angle of attack (0-5)
%   6. Mach number (0.1-0.6)

    % making the latin hypercube sample
    sample_matrix = lhsdesign(samples, numParams);
    
    % making the sampled airfoil names (NACA)
    sample_matrix(:,1) = 6*sample_matrix(:,1);  %samples for max camber, 0-6%
    sample_matrix(:,2) = 4*(sample_matrix(:,2)) + 2;   %samples for loc of max camber, 20-60% (shown as 2-6)
    
    if numParams > 2
        sample_matrix(:,3) = 6*(sample_matrix(:,3)) + 6;    %samples for max thickness, 6-12%
    end
    
    if numParams > 3
        sample_matrix(:,4) = 2*(sample_matrix(:,4)) + 5; %reynolds #, 5-7 (later converted to 1e5-1e7)
    end
    
    if numParams > 4
        sample_matrix(:,5) = 5*sample_matrix(:,5);   %alpha, 0-5degrees
    end
    
    if numParams > 5
%         sample_matrix(:,6) = ((5*sample_matrix(:,6)+1)/10).^2;    %mach, 0.1-0.6, squared
        sample_matrix(:,6) = (5*sample_matrix(:,6)+1)/10;
    end

end


        