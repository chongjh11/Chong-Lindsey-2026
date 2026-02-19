%% Make laplacian smoothing
% [L] = laplacian1D_JHC(n)
%   Note that the 1,1 and end,end has a value of "-2", you can change it to
%   zeros outside of this script
%
% - n is the number of points you want (size of matrix)
% - L is the output of the smoothing with n x x
%
%   Last modified on 3-Jul-2025
%   by Jeng Hann, Chong


function [L] = laplacian1D_JHC(n)

% n = 811;                      % Number of points
m = n+2; % adding one to expand the matrix then we will reduce to the correct size after
L = zeros(m, m);              % Initialize Laplacian matrix

% Fill in the second derivative 
for i = 2:m-1
    L(i, i-1) = 1;
    L(i, i)   = -2;
    L(i, i+1) = 1;
end

L = L(2:end-1,2:end-1); % reducing to the initial input size;

% Optional: apply zero rows to boundaries (Neumann or Dirichlet BCs)
% For Dirichlet (fixed), leave first and last row as all zeros.
% For Neumann (free slope), you could try this:
% L(1,1:2) = [-1, 1];      % Forward difference at first point
% L(n,n-1:n) = [1, -1];    % Backwards

% Applying the first and final row correctly with the 2nd and end-1 rows
L(end,end-2:end)=L(end-1,end-2:end);
L(1,1:4)=L(2,1:4);

end
