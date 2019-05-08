function [Stats, Raw, Lim] = readres(filename);
    
% CT/RJ-Mix Transdimensional MCMC for Gaussian mixtures
%
% readres: Reads the results of rj_mix in OCTAVE/MATLAB
%
% Stats = readres('filename') or
%   [Stats, Raw, Lim] = readres('filename')
% where
%   Raw(Lim(n,1):Lim(n,2)) contains w(1), mu(1), var(1),
%   ..., w(k(n)), mu(k(n)), var(k(n)) for iteration n.
%
% $Id: readres.m,v 1.1 2003/11/28 16:50:19 cappe Exp $
%
% Copyright (C) 2003, Olivier Cappé, Tobias Rydén, Christian P. Robert

eval(['Stats = load(''' filename '.sts'');']);

if (nargout > 1)
  % Raw data
  [fid,message] = fopen([filename '.res'], 'r');
  if (fid == -1)
    error(message);
  else
    Raw = fread(fid, Inf, 'double');
    fclose(fid);
    ndat = length(Raw)/3;
    if (ndat ~= floor(ndat))
      error('Doesn''t seem to be a correct result file\n');
    end
  end
  % Limits
  L = cumsum(Stats(:,2));
  nit = size(Stats,1);
  Lim = [[1; (1+3*L(1:nit-1))] , 3*L];
  if (Lim(nit,2) ~= length(Raw))
    error('Stat file and result file don''t match\n');
  end
end
