Npoints = 5000; %number of points in the metastable states
baseint = 100; %intensity of state 1
sigmas = [0 5 20]; %variance of the metastable states
Nsegmentseach = 10; %number of segments in each metastable state
traj = [];
cltraj = [];
%rng(0);

for sigma=sigmas
for percentminor=0.50:0.05:0.50 %change these numbers to get different population fraction of state 2
    for intratio=2.0:0.05:2.0 %change these numbers to get different "distance" between states
        mu2 = baseint*intratio;
        for tl=1:1:100 %change these numbers to get different transition lengths
        traj = [];
        cltraj = [];
        totalminor = percentminor*Npoints;
        totalmajor = Npoints-totalminor;
        %randomly generate segment lengths
        
        lengthsminor = round(randfixedsum(Nsegmentseach,1,totalminor,100,totalminor-100*Nsegmentseach));
        lengthsmajor = round(randfixedsum(Nsegmentseach,1,totalmajor,100,totalmajor-100*Nsegmentseach));
        if sum(lengthsminor)+sum(lengthsmajor)>Npoints
          lengthsminor(Nsegmentseach) = lengthsminor(Nsegmentseach)-(sum(lengthsminor)+sum(lengthsmajor)-Npoints);
        elseif sum(lengthsminor)+sum(lengthsmajor)<Npoints
          lengthsminor(Nsegmentseach) = lengthsminor(Nsegmentseach)+(Npoints-sum(lengthsminor)-sum(lengthsmajor));
        end
        
        
        %if you want equal segment lengths, uncomment the lines below and
        %change the 50 to the appropriate number of transition segments
        % lengthsminor = ones(50,1)*100;
        % lengthsmajor = ones(50,1)*100;
        for i=1:Nsegmentseach
            traj = vertcat(traj,laprnd(lengthsmajor(i),1,baseint,sigma));
            changes(4*i-3,1) = length(traj); 
            cltraj = vertcat(cltraj,ones(lengthsmajor(i),1));
            tlx = tl;
            for j = 1:tlx
                lambda = (tlx+1-j)/(tlx+1);
                traj = vertcat(traj,laprnd(1,1,(1-lambda)*(mu2-baseint)+baseint,sigma));
                cltraj = vertcat(cltraj,ones(1,1)*2);
            end
            changes(4*i-2,1) = length(traj);
            traj = vertcat(traj,laprnd(lengthsminor(i),1,baseint*intratio,sigma));
            cltraj = vertcat(cltraj,ones(lengthsminor(i),1)*1);
            changes(4*i-1,1) = length(traj);
            if i<Nsegmentseach
            tlx = tl;
            for j = 1:tl
                lambda = (tlx+1-j)/(tlx+1);
                traj = vertcat(traj,laprnd(1,1,(1-lambda)*(baseint-mu2)+mu2,sigma));
                cltraj = vertcat(cltraj,ones(1,1)*2);
            end
            end
            changes(4*i,1) = length(traj);
        end
                
        writematrix(traj,strcat('sigma_',num2str(sigma),'_transition_',...
            num2str(tl),'.txt'),'Delimiter',' ');
        writematrix(vertcat([0],changes(1:length(changes)-1)),strcat('sigma_',num2str(sigma),'_transition_',...
            num2str(tl),'_GroundTruth.txt'),'Delimiter',' ');
       
         gscatter(1:length(traj),traj,cltraj,colormap(jet(2)));
         ylim([min(traj)-50 max(traj)+50]);
         xlim([0 length(traj)]);
         ax = gca;
         ax.FontSize = 24;
         ax.FontWeight = 'bold';
         set(gcf,'Position',[0 0 1600 300]);
         xlabel('Time step');
         ylabel('Intensity');
         legend('off');
         for i = 1:length(changes)
             xline(changes(i));
         end
         %saveas(gcf,'traj.png');
         close(gcf);
        end
    end
end
end

function y  = laprnd(m, n, mu, sigma)
%LAPRND generate i.i.d. laplacian random number drawn from laplacian distribution
%   with mean mu and standard deviation sigma. 
%   mu      : mean
%   sigma   : standard deviation
%   [m, n]  : the dimension of y.
%   Default mu = 0, sigma = 1. 
%   For more information, refer to
%   http://en.wikipedia.org./wiki/Laplace_distribution
%   Author  : Elvis Chen (bee33@sjtu.edu.cn)
%   Date    : 01/19/07
%Check inputs
if nargin < 2
    error('At least two inputs are required');
end
if nargin == 2
    mu = 0; sigma = 1;
end
if nargin == 3
    sigma = 1;
end
% Generate Laplacian noise
u = rand(m, n)-0.5;
b = sigma / sqrt(2);
y = mu - b * sign(u).* log(1- 2* abs(u));
end


function [x,v] = randfixedsum(n,m,s,a,b)
% [x,v] = randfixedsum(n,m,s,a,b)
%
%   This generates an n by m array x, each of whose m columns
% contains n random values lying in the interval [a,b], but
% subject to the condition that their sum be equal to s.  The
% scalar value s must accordingly satisfy n*a <= s <= n*b.  The
% distribution of values is uniform in the sense that it has the
% conditional probability distribution of a uniform distribution
% over the whole n-cube, given that the sum of the x's is s.
%
%   The scalar v, if requested, returns with the total
% n-1 dimensional volume (content) of the subset satisfying
% this condition.  Consequently if v, considered as a function
% of s and divided by sqrt(n), is integrated with respect to s
% from s = a to s = b, the result would necessarily be the
% n-dimensional volume of the whole cube, namely (b-a)^n.
%
%   This algorithm does no "rejecting" on the sets of x's it
% obtains.  It is designed to generate only those that satisfy all
% the above conditions and to do so with a uniform distribution.
% It accomplishes this by decomposing the space of all possible x
% sets (columns) into n-1 dimensional simplexes.  (Line segments,
% triangles, and tetrahedra, are one-, two-, and three-dimensional
% examples of simplexes, respectively.)  It makes use of three
% different sets of 'rand' variables, one to locate values
% uniformly within each type of simplex, another to randomly
% select representatives of each different type of simplex in
% proportion to their volume, and a third to perform random
% permutations to provide an even distribution of simplex choices
% among like types.  For example, with n equal to 3 and s set at,
% say, 40% of the way from a towards b, there will be 2 different
% types of simplex, in this case triangles, each with its own
% area, and 6 different versions of each from permutations, for
% a total of 12 triangles, and these all fit together to form a
% particular planar non-regular hexagon in 3 dimensions, with v
% returned set equal to the hexagon's area.
%
% Roger Stafford - Jan. 19, 2006
% Check the arguments.
if (m~=round(m))|(n~=round(n))|(m<0)|(n<1)
 error('n must be a whole number and m a non-negative integer.')
elseif (s<n*a)|(s>n*b)|(a>=b)
 error('Inequalities n*a <= s <= n*b and a < b must hold.')
end
% Rescale to a unit cube: 0 <= x(i) <= 1
s = (s-n*a)/(b-a);
% Construct the transition probability table, t.
% t(i,j) will be utilized only in the region where j <= i + 1.
k = max(min(floor(s),n-1),0); % Must have 0 <= k <= n-1
s = max(min(s,k+1),k); % Must have k <= s <= k+1
s1 = s - [k:-1:k-n+1]; % s1 & s2 will never be negative
s2 = [k+n:-1:k+1] - s;
w = zeros(n,n+1); w(1,2) = realmax; % Scale for full 'double' range
t = zeros(n-1,n);
tiny = 2^(-1074); % The smallest positive matlab 'double' no.
for i = 2:n
 tmp1 = w(i-1,2:i+1).*s1(1:i)/i;
 tmp2 = w(i-1,1:i).*s2(n-i+1:n)/i;
 w(i,2:i+1) = tmp1 + tmp2;
 tmp3 = w(i,2:i+1) + tiny; % In case tmp1 & tmp2 are both 0,
 tmp4 = (s2(n-i+1:n) > s1(1:i)); % then t is 0 on left & 1 on right
 t(i-1,1:i) = (tmp2./tmp3).*tmp4 + (1-tmp1./tmp3).*(~tmp4);
end
% Derive the polytope volume v from the appropriate
% element in the bottom row of w.
v = n^(3/2)*(w(n,k+2)/realmax)*(b-a)^(n-1);
% Now compute the matrix x.
x = zeros(n,m);
if m == 0, return, end % If m is zero, quit with x = []
rt = rand(n-1,m); % For random selection of simplex type
rs = rand(n-1,m); % For random location within a simplex
s = repmat(s,1,m);
j = repmat(k+1,1,m); % For indexing in the t table
sm = zeros(1,m); pr = ones(1,m); % Start with sum zero & product 1
for i = n-1:-1:1  % Work backwards in the t table
 e = (rt(n-i,:)<=t(i,j)); % Use rt to choose a transition
 sx = rs(n-i,:).^(1/i); % Use rs to compute next simplex coord.
 sm = sm + (1-sx).*pr.*s/(i+1); % Update sum
 pr = sx.*pr; % Update product
 x(n-i,:) = sm + pr.*e; % Calculate x using simplex coords.
 s = s - e; j = j - e; % Transition adjustment
end
x(n,:) = sm + pr.*s; % Compute the last x
% Randomly permute the order in the columns of x and rescale.
rp = rand(n,m); % Use rp to carry out a matrix 'randperm'
[ig,p] = sort(rp); % The values placed in ig are ignored
x = (b-a)*x(p+repmat([0:n:n*(m-1)],n,1))+a; % Permute & rescale x
return
end