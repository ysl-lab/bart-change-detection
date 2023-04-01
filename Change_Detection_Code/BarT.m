%INPUT: Variable named 'tl' - vector, type: double

InitGuessFactors = [25];
maxIter = 20;
rng(0);
%SET TUNING PARAMETER VALUE HERE
lambda = 20;
chunk_to_use = length(tl);

penalties = zeros(size(tl,1),1);
for t = 1:size(tl,1)
    penalties(t) = lambda*random('Uniform',0.9,1);
end

 shift_and_merge = false;
 prev_changes = [];
 changes = find_changes(tl(1:chunk_to_use), InitGuessFactors, lambda, penalties,...
     length(tl(1:chunk_to_use)), [0]);


 changes = changes{size(tl,1),1};
 changes = table2array(cell2table(changes));
 changes(length(changes)+1) = length(tl);
 changes = transpose(changes);

iters = 0;
while true
    if shift_and_merge
        i = 2;
        while i < length(changes)
            disp(i);
            maxLike = -Inf;
            maxind = changes(i);
            for j = changes(i-1):changes(i+1)
                dummy_changes = changes;
                dummy_changes(i) = j;
                merged = false;
                if(j==changes(i-1) || j==changes(i+1))
                    dummy_changes(i) = [];
                    merged = true;
                end

                %2 neighbors + two segments governed by shifting CP
                likeN1 = 0;
                likeN2 = 0;
                likeS1 = 0;
                likeS2 = 0;

                if i > 2
                    if i == 3
                        likeN1 = find_likelihood(tl(dummy_changes(1)+1:dummy_changes(2)));
                    else
                        [n1m,n1t] = find_likelihood(tl(dummy_changes(i-2)+1:dummy_changes(i-1)),...
                        tl(dummy_changes(i-3)+1:dummy_changes(i-2)),...
                        tl(dummy_changes(i-1)+1:dummy_changes(i)));
                    likeN1 = max(n1m,n1t);
                    end
                end

                if i < length(changes)-1
                    if i == length(changes)-2
                        if merged == false
                            likeN2 = find_likelihood(tl(dummy_changes(i+1)+1:dummy_changes(i+2)));
                        else
                            likeN2 = find_likelihood(tl(dummy_changes(i)+1:dummy_changes(i+1)));
                        end
                    else
                        if merged == false
                            [n2m,n2t] = find_likelihood(tl(dummy_changes(i+1)+1:dummy_changes(i+2)),...
                                tl(dummy_changes(i)+1:dummy_changes(i+1)),...
                                tl(dummy_changes(i+2)+1:dummy_changes(i+3)));
                            likeN2 = max(n2m,n2t);
                        else
                            [n2m,n2t] = find_likelihood(tl(dummy_changes(i)+1:dummy_changes(i+1)),...
                                tl(dummy_changes(i-1)+1:dummy_changes(i)),...
                                tl(dummy_changes(i+1)+1:dummy_changes(i+2)));
                            likeN2 = max(n2m,n2t);
                        end
                    end
                end

                if i > 2 && (i < length(dummy_changes) || merged == false)
                    [s1m,s1t] = find_likelihood(tl(dummy_changes(i-1)+1:dummy_changes(i)),...
                        tl(dummy_changes(i-2)+1:dummy_changes(i-1)),...
                        tl(dummy_changes(i)+1:dummy_changes(i+1)));
                    likeS1 = max(s1m,s1t);
                elseif i <= 2
                    likeS1 = find_likelihood(tl(dummy_changes(i-1)+1:dummy_changes(i)));
                else
                    likeS1 = find_likelihood(tl(dummy_changes(i-1)+1:dummy_changes(i)));
                end

                if i < length(changes)-1 && merged == false
                    [s2m,s2t] = find_likelihood(tl(dummy_changes(i)+1:dummy_changes(i+1)),...
                        tl(dummy_changes(i-1)+1:dummy_changes(i)),...
                        tl(dummy_changes(i+1)+1:dummy_changes(i+2)));
                    likeS2 = max(s2m,s2t);
                elseif i < length(changes) && merged == false
                    likeS2 = find_likelihood(tl(dummy_changes(i)+1:dummy_changes(i+1)));
                end

                currentLike = likeN1 + likeN2 + likeS1 + likeS2 - lambda*length(dummy_changes);
                if currentLike > maxLike
                    maxLike = currentLike;
                    maxind = j;
                end
            end
            changes(i) = maxind;
            if(maxind==changes(i-1) || maxind==changes(i+1))
                    changes(i) = [];
                    i = i-1;
            end
            i = i+1;
        end
    end

    if isequal(prev_changes,changes) || iters == maxIter
        disp("HELLO")
        break;
    end

    prev_changes = changes;
    if length(changes) == length(prev_changes)
        shift_and_merge = true;
    end

    iters = iters + 1;
end


[likesum,tbl] = sum_like(tl,InitGuessFactors,changes,lambda);
writematrix(changes,'changes.txt','Delimiter',' ');
grouping = zeros(length(tl),1);
for i = 1:length(changes)-1
    if(tbl(i,1)==1)
        grouping(changes(i)+1:changes(i+1)) = 1;
    elseif(tbl(i,1)==2)
        grouping(changes(i)+1:changes(i+1)) = 2;
    else
        grouping(changes(i)+1:changes(i+1)) = 3;
    end
end
 scatter(find(grouping==1),tl(find(grouping==1)),'.r');
 hold on;
 scatter(find(grouping==2),tl(find(grouping==2)),'.b');
 for i = 2:length(changes)
 xline(changes(i));
 end
%  limits =  ylim;
%  for i = 2:length(gtchanges)-1
%      scatter(gtchanges(i),limits(1),24,'o','MarkerEdgeColor','k','MarkerFaceColor','y');
%  end
ax = gca;
ax.FontSize = 24;
ax.FontWeight = 'bold';
xlabel('Time step')
ylabel('Value')
set(gcf,'Position',[0 0 800 200]);
xlim([0 length(tl)]);

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

function [sumlike,tbl] = sum_like(tl,~,changelist,lambda)
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
    sumlike = sumlike + max(meta_seg,trans_seg);
    if meta_seg > trans_seg
        tbl(i,1) = 1;
    else
        tbl(i,1) = 2;
    end
    tbl(i,2) = meta_seg;
    tbl(i,3) = trans_seg;
end

sumlike = sumlike-lambda*(length(changelist)-2);
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