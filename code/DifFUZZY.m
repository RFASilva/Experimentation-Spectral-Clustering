function U = DifFUZZY(dat, M, gammas)

% DifFUZZY clusters a dataset "dat" using the algorithm presented in
% "DifFUZZY: A fuzzy spectral clustering algorithm for complex data sets"
%
% inputs:
%     dat - N*p matrix, input data (N: number of data points, p: dimension
%           of the data points)
%     M - positive integer <= N, required parameter related to the minimum
%         core clusters size
%     gammas - 3 vector, optional parameters, gamma_1, gamma_2 and gamma_3
%              Intuitively,
%              gammas(1): time scale of random walk          (default: 0.3)
%              gammas(2): length of time step of random walk (default: 0.1)
%              gammas(3): number of steps of random walk     (default: 1.0)
% output:
%
%     U - N*C Membership values matrix, where C: number of clusters, found
%         by DifFUZZY
%
% Example
%       noise = 0.2;
%       x1 = ones(40,1) + 1*noise*randn(40,1);
%       y1 = ones(40,1) + 1*noise*randn(40,1);
%       x2 = 2*ones(120,1) + 1.5*noise*randn(120,1);
%       y2 = 2.1*ones(120,1) + 1.5*noise*randn(120,1);
%       x = [x1; x2];
%       y = [y1; y2];
%       dat = [x y];
%       U = DifFUZZY(dat, 30, [0.3; 0.1; 1.0]);
%
% Copyright (c) 2009 O. Cominetti, A. Matzavinos, D. Kulasiri,
%                    S. Samarasinghe, S. Liu, P.K. Maini and R. Erban
%
% DifFUZZY toolbox version 1.1
% This code is distributed under the terms of the GNU General Public
% License.  See "license.txt" for details.

if nargin ~= 2 & nargin ~= 3,
  display('Too few or too many input arguments');
  return
end

gammasDefault = [0.3, 0.1, 1];

if nargin == 2
  gamma_1 = gammasDefault(1);
  gamma_2 = gammasDefault(2);
  gamma_3 = gammasDefault(3);
end

if nargin == 3
  
  for i = 1:3
    if isnan(gammas(i))
      gammas(i) = gammasDefault(i);
    end
  end
  
  gamma_1 = gammas(1);
  gamma_2 = gammas(2);
  gamma_3 = gammas(3);
  
end

m = length(dat(:,1));
NumDim = length(dat(1,:));

%--------------------------------------------------------------------------

[sid, NC, rr, nclsgt, clusts] = FindSigmaStar(dat, M);

if NC == 0
  display('No core clusters');
  return
end

if sum(clusts == 0) == 0
  undecIdx = []; 
  U = [];
  display('No undecided points');
else
  
  undecIdx = zeros(sum(clusts==0),1);
  g = 1;
  for mm = 1:length(clusts)
    if clusts(mm) == 0
      undecIdx(g) = mm;
      g = g+1;
    end
  end
  
end

%--------------------------------------------------------------------------

display('Finding beta^*')
display(' ')

[betaStar, betas, L] = FindBetaStar(dat, clusts, NC, gamma_1);

%--------------------------------------------------------------------------

U = zeros(length(clusts), NC);

for cl = 1:length(clusts)
  for nc = 1:NC
    U(cl,nc) = (clusts(cl) == nc);
  end
end

display(['Number of undecided points:', num2str(length(undecIdx))]);
display(' ');
display(['Number of clusters:', num2str(NC)]);
display(' ');

%--------------------------------------------------------------------------

display('Generating the weight matrix')
display(' ')

clear W;
W = zeros(m,m);

W = exp((-(distance(dat',dat')).^2)./betaStar);

for i = 1:m
  for j = 1:m
    if ((clusts(i)~=0) & (clusts(i)~=0))
      if ((clusts(i)==clusts(j)) & (clusts(i)~=0))
        W(i,j) = 1;
      end
    end
  end
end

%--------------------------------------------------------------------------

display('Evaluating membership functions')
display(' ')

D = diag(sum(W));
P = eye(m,m) + (W - D).*gamma_2/max(max(D));
V = eigs(P,2);
lambda2 = real(V(2));
alpha = floor(gamma_3/abs(log(lambda2)));

if sum(sum((P^alpha) <= 0)) ~= 0
  display('Error -> Negative entries in transition matrix P^alpha')
  return
end


for i = 1:length(undecIdx)
  s=undecIdx(i);
  memf = FindDiffDist(dat, W, s, clusts, NC, P, gamma_2, alpha);
  U(s,:) = memf;
end

%--------------------------------------------------------------------------

display('Generating HCT plots')
display(' ')

HCTPlot(clusts,U,dat,undecIdx,NC,0.9);

%--------------------------------------------------------------------------

display('Generating Membership Values plot')
display(' ')

MVPlot(U);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sid,NC,rr,nclsgt,nclust]= FindSigmaStar(dat, M)

a = max(max(distance(dat',dat')));
dat = dat./a;
sigmaEnd = 1;

m = length(dat(:,1));
ncls = [];
ii = 0;

for sigma = 0.001:0.005:sigmaEnd
  ii=ii+1;
  
  clear W;
  W = zeros(m,m); % initialize weight matrix
  W = ((distance(dat',dat')) < sigma);
  W = W - diag(ones(m,1));
  
  vcii = zeros(m,1);
  
  z = cmp(W,m);
  
  clear zz
  zz = z(:,1);
  
  for kk=2:1:size(z,2)
    zz = zz + kk*z(:,kk);
  end
  
  vcii(1:length(sum(z,1)),1) = sum(z,1)';
  Vv(:,ii) = vcii;
  
end

clear Vv2;

for i = 1:1:size(Vv,2)
  for j = 1:1:m
    nm = find(Vv(:,i) == j);
    if (length(nm)>0)
      Vv2(j,i) = length(nm);
    end
  end
end

nclsgt = [];

for i = 1:1:size(Vv2,2)
  nclsgt = [nclsgt sum(Vv2(M:end,i))];
end

dncl = find(nclsgt == max(nclsgt));

dncl(1);

rr = 0.001:0.005:sigmaEnd;

sid = rr(dncl(1));
NC = max(nclsgt);

% here we determine to which cluster each element belongs to.
W2 = ((distance(dat',dat')) < sid);

W2 = W2 - diag(ones(m,1));
z = cmp(W2,m);

nclust = z(:,1);

for kk = 2:1:size(z,2)
  nclust = nclust + kk*z(:,kk);
end

% now we put a zero in nclust for each point that does not belong to a 
% cluster (that way we will be able to identify it as an undecided point)
for i = 1:max(nclust)
  if sum(nclust==i) < M
    nclust = nclust.*(nclust~=i);
  end
end

% now we make sure the cluster numbers begin from 1 and are on assending
% order.
vsdz = zeros(NC,1);
k = 1;

% first we identify the cluster numbers
for i = 1:length(nclust)
  if (nclust(i) ~= 0)
    if ((sum(vsdz == nclust(i)) == 0) == 1)
      vsdz(k) = nclust(i);
      k = k+1;
    end
  end
end


for i=1:length(vsdz)
  nclust = nclust.*(nclust~=vsdz(i)) + i.*(nclust==vsdz(i)).*(nclust~=0);
end

sid = sid*a;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [betaStar,betas,L,a] = FindBetaStar(dat, clusts, NC, gamma_1)

a = max(max(distance(dat',dat')));
dat = dat./a;

betas = logspace(-7,4);
betanum = length(betas);
m = length(dat(:,1));
DistMat = -1*(distance(dat',dat').^2);

L = zeros(betanum,1);

for beta=1:betanum
  
  clear W;
  W = zeros(m,m);
  W = exp(DistMat./betas(beta));
  
  B = repmat(clusts,[1 length(clusts)]);
  A = B==B';
  W(logical((B~=0).*A)) = 1;
  
  L(beta) = sum(sum(W));
  
end

LStar = min(L) + gamma_1*(max(L) - min(L));

for h = 1:length(L)
  if LStar < L(h)
    break;
  else
  end
end

betaStar=(a^2).*betas(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function memf= FindDiffDist(dat, W, s, clusts, NC, P, gamma_2, alpha)

m = length(dat(:,1));
pclos_idx(1:NC) = s;
distmin(1:NC) = realmax;

D = diag(sum(W));
DeltaT = 0.1/max(max(D));


for nc = 1:NC
  for k = 1:m
    if clusts(k) == nc
      if distance(dat(s,:)',dat(k,:)') < distmin(nc)
        pclos_idx(nc) = k;
        distmin(nc) = distance(dat(s,:)',dat(k,:)');
      end
    end
  end
end

uu = W(s,:);
vv = W(:,s);

e = zeros(m,1);
e(s) = 1;

PP = P^alpha;
PP = PP*e;

for nc = 1:NC
  
  W(s,:) = W(pclos_idx(nc),:);
  W(:,s) = W(:,pclos_idx(nc));
  W(s,s) = 1;
  
  D = diag(sum(W));
  D = W-D;
  D = D.*DeltaT;
  
  
  P2 = eye(m,m) + D;
  
  W(s,:) = uu;
  W(:,s) = vv;
  
  clear D;
  
  PP2 = P2^alpha;
  PP2 = PP2*e;
  
  DL0 = sum((PP-PP2).^2);
  DL = sqrt(DL0);
  clear DL0;
  
  Dist(nc) = DL;
  
  clear DL;
  
end

memf = (1./Dist)/sum(Dist.^(-1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HCTPlot(clusts, U, dat, undecIdx, NC, z)

NumDim = length(dat(1,:));
clusts2 = clusts;

for nc = 1:NC
  for i = 1:length(undecIdx)
    if U(undecIdx(i),nc) > z
      clusts2(undecIdx(i)) = nc;
    end
  end
end

n = length(U(1,:));

if NC==1
    if isempty(undecIdx)
        kmap = [1,0,0.9;];
    else
        kmap = [0,0,0; 1,0,0.9;];

    end
end

if NC==2
    if isempty(undecIdx)
        kmap = [1,0,0.9; 0,1,0; 0,0.5,1];
    else
        kmap = [0,0,0; 1,0,0.9; 0,1,0; 0,0.5,1];
    end
end

if NC==3
    if isempty(undecIdx)
        kmap = [0,1,0; 1,0,0.9; 0,0.5,1];
    else
        kmap = [0,0,0; 0,1,0; 1,0,0.9; 0,0.5,1];

    end
end

if NC==4
    if isempty(undecIdx)
        kmap = [0,1,0; 0,1,0; 1,0,0.9; 0,0.5,1];
    else
        kmap = [0,0,0; 0,1,0; 1,0,0.9; 0,1,1; 0,0.5,1];

    end
end    

if NC==5
    if isempty(undecIdx)
        kmap = [0,1,0; 1,1,0; 1,0,0.9; 0,1,1; 0,0.5,1];
    else
        kmap = [0,0,0; 0,1,0; 1,1,0; 1,0,0.9; 0,1,1; 0,0.5,1];
    end
end  

if NumDim == 2
  figure()
  colormap(kmap);
  scatter(dat(:,1),dat(:,2),18,clusts2,'filled');
  xlabel('x','FontSize',14)
  ylabel('y','FontSize',14)
  axis tight
  axis equal
  h=gca;
  set(h,'FontSize',14);
  set(h,'box','on');
  title('HCT Plot');
end

if NumDim == 3
  figure()
  colormap(kmap);
  scatter3(dat(:,1),dat(:,2),dat(:,3),18,clusts2,'filled');
  xlabel('x','FontSize',14);
  ylabel('y','FontSize',14);
  zlabel('z','FontSize',14);
  axis tight
  axis equal
  h=gca;
  set(h,'FontSize',14);
  set(h,'box','on');
  title('HCT Plot');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MVPlot(U)

m = length(U(:,1));

kmap = [0,1,0; 1,0,0.9; 0,0.5,1; 1,1,0; 0,1,1];

figure()
bar(U,'stacked'), colormap(kmap);
axis([0.5 (m+0.5) 0 1]);
ylabel('Membership values','FontSize',14);
xlabel('Data Point Number','FontSize',14);
h=gca;
set(h,'FontSize',14);
title('MV Plot') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = distance(a,b)

% DISTANCE - computes Euclidean distance matrix
% Author: Roland Bunschoten, University of Amsterdam

if (size(a,1) ~= size(b,1))
  error('A and B should be of same dimensionality');
end

aa = sum(a.*a,1); 
bb = sum(b.*b,1); 
ab = a'*b;
d = sqrt(abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = cmp(W,m)

z = false(m,1);
nnn = true(1,m);

while any(nnn)
  cm = false(1, m); % marks nodes in the same component as node 1
  N = false(1, m);
  k = min(find(nnn == true));
  N(k) = true;
  while any(N)
    cm = cm | N;
    nnn = nnn & (~N);
    N = sum(W(N,:), 1) & ~cm;
  end
  z = [z cm'];
end

z=z(:,2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




