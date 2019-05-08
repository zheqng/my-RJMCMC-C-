function [Stats, Raw, Lim] = readres2(filename);

% CT/RJ-Mix Transdimensional MCMC for Gaussian mixtures
%
% readres2: Reads the results of ct_mix in OCTAVE/MATLAB
%
% Stats = readres2('filename') or
%   [Stats, Raw, Lim] = readres2('filename')
% where
%   Raw(Lim(n,1):Lim(n,2)) contains wght, w(1), mu(1), var(1),
%   ..., w(k(n)), mu(k(n)), var(k(n)) for iteration n.
%
% $Id: readres2.m,v 1.1 2003/11/28 16:50:19 cappe Exp $
%
% Copyright (C) 2003, Olivier Cappé, Tobias Rydén, Christian P. Robert

eval(['Stats = load(''' filename '.st2'');']);

if (nargout > 2)
  % Raw data
  [fid,message] = fopen([filename '.rs2'], 'r');
  if (fid == -1)
    error(message);
  else
    Raw = fread(fid, Inf, 'double');
    fclose(fid);
  end
  % Limits
  L = cumsum(Stats(:,2));
  nit = size(Stats,1);
  G = 3*L + (1:nit)';
  Lim = [[1; (1+G(1:nit-1))] , G];
  if (Lim(nit,2) ~= length(Raw))
    error('Stat file and result file don''t match\n');
  end
end
