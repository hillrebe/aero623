function avg = roeavg(rhoL, rhoR, varL, varR)

num = sqrt(rhoL)*varL + sqrt(rhoR)*varR;
den = sqrt(rhoL) + sqrt(rhoR);
avg = num/den;

end