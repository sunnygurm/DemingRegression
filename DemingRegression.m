% Written by Sunny Gurm (c) 2018- credit to James Hall; his code was used to calculate
% confidence intervals based on the jacknife method. Preamble from his code is 
% below, as per his conditions:

% Copyright (c) 2011, Hidden Solutions, LLC
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%   * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
%    * Neither the name Hidden Solutions, LLC nor the names of any
%      contributors may be used to endorse or promote products derived from 
%      this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL HIDDEN SOLUTIONS, LLC BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%


function [beta, stats] = DemingRegression(x,y,lambda,alpha)

  n = length(x);
  xbar = mean(x);
  ybar = mean(y);
  covmat = cov([x',y'])
  Sxx = (1/(n-1)) * sum((x - xbar).^2)
  Sxy = (1/(n-1)) *sum((x - xbar).*(y-ybar));
  Syy = (1/(n-1)) * sum((y - ybar).^2)
  
  B1 = (Syy - lambda*Sxx + sqrt((Syy-lambda*Sxx)^2 + 4*lambda*Sxy^2))/(2*Sxy);
  
  B0 = ybar-B1*xbar;
  
  beta = [B0; B1];
  
  %compute conf intervals below
  b_sub = zeros(2,n);
  ignoreFlag = [false; true(n-1,1)];
  for nn = 1:n
    b_sub(:,nn) = Deming(x(circshift(ignoreFlag,nn)),y(circshift(ignoreFlag,nn)),lambda);
  end
    
  % Compute standard error of the slope and intercept
  stats.s_b = std(b_sub,[],2)*(n-1)/sqrt(n);
  
  % Statistics toolbox is installed
  stats.t_c    = tinv(1-alpha/2,n-2);
  stats.b_ci   = beta*[1 1] + stats.t_c*stats.s_b*[-1 1];
  
  
end