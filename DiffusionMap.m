function [mPhi, mLam] = DiffusionMap(mY)
% DiffusionMap - Create Diffusion Map from Data mY
%   mY is a matrix, whos rows are individual observations (e.g. frames),
%   and its columns are individual observed variables (e.g. pixels).
mW         = squareform( pdist(mY) );
eps        = median(mW(:));
mK         = exp(-mW.^2 / (eps^2));
mA         = mK ./ sum(mK, 2);

[mPhi, mLam] = eig(mA);
end