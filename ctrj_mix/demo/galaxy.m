% CT/RJ-Mix Transdimensional MCMC for Gaussian mixtures
%
% demo.m: Example of script that reads the results of rj_mix and ct_mix in
% MATLAB to plot the estimated posterior for k (number of components) and
% density estimate (note: the script works in OCTAVE as well except for the
% plots).
%
% $Id: galaxy.m,v 1.1 2003/11/28 16:50:19 cappe Exp $
%
% Copyright (C) 2003, Olivier Cappé, Tobias Rydén, Christian P. Robert

% The function gauseval is taken from the H2M toolbox
% (http://www.tsi.enst.fr/~cappe/h2m)
path(path, 'h2m');

% Load the data
X = load('galaxy.dat');
NGRID = 50;
t = linspace(-1.3,1.3,50)';
M = 15;

%%% RJ Sampler
[Stats, Raw, Lim] = readres('galaxy');
DHIST = zeros(M,1);
DENS = zeros(NGRID,1);
for i = 2:size(Stats,1)       % Omit initial value
  % Histogram for k
  DHIST(Stats(i,2)) = DHIST(Stats(i,2)) + 1;  
  % Density estimate
  R = reshape(Raw(Lim(i,1):Lim(i,2)), Stats(i,2), 3);
  w = R(:,1);
  mu = R(:,2);
  var = R(:,3);
  dens = gauseval(t, mu, var, 1);
  DENS = DENS + dens*w;
end
DHIST = DHIST/(size(Stats,1)-1);
DENS = DENS/(size(Stats,1)-1);
% Plots
clf;
subplot(221);
bar(1:10, DHIST(1:10), 1);
axis([0 10.5 0 0.3]);
ylabel('reversible jump sampler');
subplot(222);
[nn,xx] = hist(X,50);
bar(xx,(nn/length(X))/(xx(2)-xx(1)), 0.01);
hold on;
plot(t, DENS, 'r');

%%% CT Sampler
[Stats2, Raw2, Lim2] = readres2('galaxy');
DHIST2 = zeros(M,1);
DENS2 = zeros(NGRID,1);
nrm = 0;
for i = 2:size(Stats2,1)       % Omit initial value
  % Weight (associated with CT duration of the simulated state)
  weight = Raw2(Lim2(i,1));
  % Histogram for k
  DHIST2(Stats2(i,2)) = DHIST2(Stats2(i,2)) + weight;
  nrm = nrm + weight;
  % Density estimate
  R = reshape(Raw2(Lim2(i,1)+1:Lim2(i,2)), Stats2(i,2), 3);
  w = R(:,1);
  mu = R(:,2);
  var = R(:,3);
  dens = gauseval(t, mu, var, 1);
  DENS2 = DENS2 + (dens*w)*weight;
end
DHIST2 = DHIST2/nrm;
DENS2 = DENS2/nrm;
% Plots
subplot(223);
bar(1:10, DHIST2(1:10), 1);
axis([0 10.5 0 0.3]);
xlabel('k');
ylabel('continuous time sampler');
subplot(224);
[nn,xx] = hist(X,50);
bar(xx,(nn/length(X))/(xx(2)-xx(1)), 0.01);
hold on;
plot(t, DENS2, 'r');
xlabel('density');
