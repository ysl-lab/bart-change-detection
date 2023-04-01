
%INPUT: Variable named 'tl' - matrix (rows are points, columns are dimensions), type: double
 InitGuessFactors = [25];
 maxIter = 20;
 rng(0);
%SET TUNING PARAMETER VALUES HERE
lambda = 20;
alpha = 0.7;
chunk_to_use = size(tl,1);

penalties = zeros(size(tl,1),1);
for t = 1:size(tl,1)
    penalties(t) = lambda*(random('Uniform',0.9,1));
end

changes_out = zeros(size(tl,1),size(tl,2));

for dim = 1:size(tl,2)
 changes_dim = find_changes(tl(1:chunk_to_use,dim), InitGuessFactors, lambda, penalties,...
     length(tl(1:chunk_to_use,:)), [0]);
 changes_dim = changes_dim{size(tl,1),1};
 changes_dim = table2array(cell2table(changes_dim));
 changes_dim(length(changes_dim)+1) = length(tl);
 changes_dim = transpose(changes_dim);
 changes_out(changes_dim(2:length(changes_dim)),dim) = 1;
end

clearvars changes_dim

for t = 1:size(tl,1)
    if ~(nnz(changes_out(t,:))==0)
        penalties(t) = lambda*(nnz(changes_out(t,:)))^alpha;
    else
        penalties(t) = 0;
    end
end
totalPen = sum(penalties);

for dim = 1:size(tl,2)
 shift_and_merge = false;
 changes_dim = vertcat([0],find(changes_out(:,dim)>0));
 prev_changes = [];
 iters = 0;
while true
    if shift_and_merge
        i = 2;
        while i < size(changes_dim,1)
            disp(i);
            maxLike = -Inf;
            maxind = changes_dim(i);
            for j = changes_dim(i-1):changes_dim(i+1)
                dummy_changes = changes_dim;
                dummy_changes(i) = j;
                merged = false;
                if(j==changes_dim(i-1) || j==changes_dim(i+1))
                    dummy_changes(i) = [];
                    merged = true;
                end
                
                likeN1 = 0;
                likeN2 = 0;
                likeS1 = 0;
                likeS2 = 0;
                
                if i > 2
                    if i == 3
                        likeN1 = find_likelihood(tl(dummy_changes(1)+1:dummy_changes(2),dim));
                    else
                        [n1m,n1t] = find_likelihood(tl(dummy_changes(i-2)+1:dummy_changes(i-1),dim),...
                        tl(dummy_changes(i-3)+1:dummy_changes(i-2),dim),...
                        tl(dummy_changes(i-1)+1:dummy_changes(i),dim));
                    likeN1 = max(n1m,n1t);
                    end
                end
                
                if i < length(changes_dim)-1
                    if i == length(changes_dim)-2
                        if merged == false
                            likeN2 = find_likelihood(tl(dummy_changes(i+1)+1:dummy_changes(i+2),dim));
                        else
                            likeN2 = find_likelihood(tl(dummy_changes(i)+1:dummy_changes(i+1),dim));
                        end
                    else
                        if merged == false
                            [n2m,n2t] = find_likelihood(tl(dummy_changes(i+1)+1:dummy_changes(i+2),dim),...
                                tl(dummy_changes(i)+1:dummy_changes(i+1),dim),...
                                tl(dummy_changes(i+2)+1:dummy_changes(i+3),dim));
                            likeN2 = max(n2m,n2t);
                        else
                            [n2m,n2t] = find_likelihood(tl(dummy_changes(i)+1:dummy_changes(i+1),dim),...
                                tl(dummy_changes(i-1)+1:dummy_changes(i),dim),...
                                tl(dummy_changes(i+1)+1:dummy_changes(i+2),dim));
                            likeN2 = max(n2m,n2t);
                        end
                    end
                end
                
                if i > 2 && (i < length(dummy_changes) || merged == false)
                    [s1m,s1t] = find_likelihood(tl(dummy_changes(i-1)+1:dummy_changes(i),dim),...
                        tl(dummy_changes(i-2)+1:dummy_changes(i-1),dim),...
                        tl(dummy_changes(i)+1:dummy_changes(i+1),dim));
                    likeS1 = max(s1m,s1t);
                elseif i <= 2
                    likeS1 = find_likelihood(tl(dummy_changes(i-1)+1:dummy_changes(i),dim));
                else
                    likeS1 = find_likelihood(tl(dummy_changes(i-1)+1:dummy_changes(i),dim));
                end
                
                if i < length(changes_dim)-1 && merged == false
                    [s2m,s2t] = find_likelihood(tl(dummy_changes(i)+1:dummy_changes(i+1),dim),...
                        tl(dummy_changes(i-1)+1:dummy_changes(i),dim),...
                        tl(dummy_changes(i+1)+1:dummy_changes(i+2),dim));
                    likeS2 = max(s2m,s2t);
                elseif i < length(changes_dim) && merged == false
                    likeS2 = find_likelihood(tl(dummy_changes(i)+1:dummy_changes(i+1),dim));
                end         
                
                if j > 0 && j < size(tl,1)
                    countConcurrent = nnz(changes_out(j,setdiff(1:size(tl,2),dim)))+1;
                    newPen = totalPen - penalties(changes_dim(i)) - ...
                        penalties(j) + lambda*(countConcurrent^alpha);
                else
                    newPen = totalPen - penalties(changes_dim(i));
                end

                currentLike = likeN1 + likeN2 + likeS1 + likeS2 - newPen;
                if currentLike > maxLike
                    maxLike = currentLike;
                    maxind = j;
                end
            end
            penalties(changes_dim(i))=0;
            if maxind > 0
            penalties(maxind)=0;
            end
            changes_dim(i) = maxind;
            if maxind > 0
            countConcurrent = nnz(changes_out(maxind,setdiff(1:size(tl,2),dim)))+1;
            penalties(maxind) = lambda*(countConcurrent^alpha);
            end
            totalPen = sum(penalties);
            if(maxind==changes_dim(i-1) || maxind==changes_dim(i+1))
                    changes_dim(i) = [];
            end
            changes_out(:,dim) = 0;
            changes_out(changes_dim(2:length(changes_dim)),dim) = 1;
            i = i+1;
        end
    end

    if isequal(prev_changes,changes_dim) || iters == maxIter
        disp("HELLO")
        break;
    end
    
    prev_changes = changes_dim;
    if length(changes_dim) == length(prev_changes)
        shift_and_merge = true;
    end
    
    iters = iters + 1;
end
end

writematrix(changes_out,"changes.txt",'Delimiter',' ');
tiledlayout(size(tl,2),1);
for dim = 1:size(tl,2)
    changes_dim = vertcat([0],find(changes_out(:,dim)>0));
    [likesum,tbl,score] = sum_like(tl(:,dim),InitGuessFactors,changes_dim,lambda);
    writematrix(tbl,strcat('log_likelihoods_',num2str(dim),'.txt'),'Delimiter',' ');
    
    %Show penalized LL
    sum_pen = 0;
    for i = 1:size(changes_out,1)
        if ~(nnz(changes_out(i,:))==0)
        sum_pen = sum_pen + lambda*(nnz(changes_out(i,:))^alpha);
        end
    end
    penLL = sum(score) - sum_pen;
    if(dim==1)
    disp(strcat("Penalized LL:",num2str(penLL)));
    disp(strcat("Penalty:",num2str(sum_pen)));
    disp(strcat("Non-penalized LL:",num2str(sum(score))));
    end
        
    
grouping = zeros(size(tl,1),1);
for i = 1:length(changes_dim)-1
    if(tbl(i,1)==1)
        grouping(changes_dim(i)+1:changes_dim(i+1)) = 1;
    elseif(tbl(i,1)==2)
        grouping(changes_dim(i)+1:changes_dim(i+1)) = 2;
    else
        grouping(changes_dim(i)+1:changes_dim(i+1)) = 3;
    end
end
nexttile;
 scatter(find(grouping==1),tl(find(grouping==1),dim),'.r');
 hold on;
 scatter(find(grouping==2),tl(find(grouping==2),dim),'.b');
 for i = 2:length(changes_dim)
 xline(changes_dim(i));
 end
ax = gca;
ax.FontSize = 24;
ax.FontWeight = 'bold';
xlabel('Time step')
ylabel('Value')
set(gcf,'Position',[0 0 800 200*size(tl,2)]);
xlim([0 length(tl)]);
end

%changes_dim1 = find(changes_out(:,1)>0);
%changes_dim2 = find(changes_out(:,2)>0);


function changes = find_changes(data,InitGuessFactors,lambda,penalties,T,s_allowed)
   F = zeros(T+1,1);
   F(1) = -lambda;
   changes{1,:} = [];
   for t = 1:T-1
       maxterm = -Inf;
       smax = 0;
       testvar = zeros(length(s_allowed),1);
       for i = 1:length(s_allowed)
           Fcurrent = F(s_allowed(i)+1);
           lmax = -Inf;
           for LSize=InitGuessFactors
           for RSize=InitGuessFactors
           [lls_meta,lls_trans] = find_likelihood(data(s_allowed(i)+1:t),...
               data(max(s_allowed(i)+1-LSize,1):s_allowed(i)),...
               data(t+1:min(t+RSize,T)));
           lls = max(lls_meta,lls_trans);
           if lls > lmax
               lmax = lls;
           end
           end
           end
           lls = lmax;
           pencurrent = penalties(s_allowed(i)+1);
           testvar(i) = Fcurrent+lls-pencurrent;
           if testvar(i) >= maxterm
               maxterm = testvar(i);
               smax = s_allowed(i);
           end
       end
       F(t+1) = maxterm;
       changes{t+1,:} = horzcat(changes{smax+1,:},{smax});
             
       s_allowed_1 = s_allowed(find(testvar >= (F(t+1) - penalties(t))));
       s_allowed_2 = s_allowed(find(testvar == -Inf));
       s_allowed = union(s_allowed_1,s_allowed_2);
       s_allowed = vertcat(s_allowed,[t]);
       
       disp(t)
   end
end

function [sumlike,tbl,score] = sum_like(tl,~,changelist,~)
segmentlist = changelist;
sumlike = 0;
tbl = zeros(size(segmentlist,1),3);
for i = 1:size(segmentlist,1)-1
    trans_seg = -Inf;
    [meta_seg, trans_seg] = find_likelihood(tl(segmentlist(i)+1:segmentlist(i+1)));
    if (i > 1 && i < size(segmentlist,1)-1)
        [meta_seg, trans_seg] = find_likelihood(tl(segmentlist(i)+1:segmentlist(i+1)),...
            tl(segmentlist(i-1)+1:segmentlist(i)),...
            tl(segmentlist(i+1)+1:segmentlist(i+2)));
    end
    score(i) = max(meta_seg,trans_seg);
    sumlike = sumlike + score(i);
    if meta_seg > trans_seg
        tbl(i,1) = 1;
    else
        tbl(i,1) = 2;
    end
    tbl(i,2) = meta_seg;
    tbl(i,3) = trans_seg;
end
end


%Main likelihood calculation function
function [metalike,translike] = find_likelihood(varargin)
epsilon = 0.000000000001;

if nargin == 3
    data = varargin{1};
    prevseg = varargin{2};
    nextseg = varargin{3};
    medianprev = median(prevseg);
    sigmaprev = std(prevseg) + epsilon;
    mediannext = median(nextseg);
    sigmanext = std(nextseg) + epsilon;
    translike = 0;
else
    data = varargin{1};
    translike = -Inf;
end

metalike = 0;
segment = data;
segmedian = median(segment);
segdev = std(segment) + epsilon;

%Compute log likelihood using general formula for both metastable and
%transition
for l = 1:length(segment)
    metalike = metalike - (log(2*segdev) + (1/segdev)...
        *abs(segment(l)-segmedian));
    if(nargin == 3)
        lambda = l/(length(segment)+1);
        mu_term = lambda*mediannext+(1-lambda)*medianprev;
        sigma = lambda*sigmanext+(1-lambda)*sigmaprev;
        translike = translike - log(2*sigma) - ((...
            abs(segment(l)-mu_term))/sigma); 
    end
end

%DIAGNOSTIC
if(length(segment)<3)
    metalike = -Inf;
    translike = -Inf;
end

end