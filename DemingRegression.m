% Written by Sunny Gurm (c) 2018

function beta = DemingRegression(x,y,lambda,alpha)

  n = length(x);
  xbar = mean(x);
  ybar = mean(y);
  covmat = cov([x',y']);
  Sxx = (1/(n-1)) * sum((x - xbar).^2);
  Sxy = (1/(n-1)) *sum((x - xbar).*(y-ybar));
  Syy = (1/(n-1)) * sum((y - ybar).^2);
  
  B1 = (Syy - lambda*Sxx + sqrt((Syy-lambda*Sxx)^2 + 4*lambda*Sxy^2))/(2*Sxy);
  
  B0 = ybar-B1*xbar;
  
  beta = [B0; B1];
  
  
end