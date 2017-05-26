function nulldist = connISC_stats(roits,iterations)
%CONNICS_COMPUTE Summary of this function goes here
%   Detailed explanation goes here

N_subj = length(roits);

for i=1:N_subj
    if abs(1-var(roits{i}(:,1)))>1e-6 || abs(mean(roits{i}(:,1)))>1e-6
        error('Data is not z-scored!')
    end        
end


N_roi = size(roits{1},2);

N_subj_pairs = N_subj^2 - N_subj;

T = size(roits{1},1);
div = T-1;

ss=sort(randsample(N_roi,iterations,true));
tt=zeros(size(ss));
uni=unique(ss);
for i=1:length(uni)
    ind=find(ss==uni(i));
    a=1:N_roi;
    a(uni(i))=[];
    tt(ind)=randsample(a,length(ind),true);
end

fprintf('Computing connectivity data statistics\n')

nulldist=zeros(1,iterations,'single');
parfor iter= 1:iterations
    
    s=ss(iter);
    t=tt(iter);
    
    if s==t
        error('Bug found!')
    end
    
    shift=randi(T,1,N_subj);
    roidata1 = zeros(T,N_subj);
    roidata2 = zeros(T,N_subj);
    for i=1:N_subj
        roidata1(:,i)=roits{i}([((shift(i)+1):end),(1:shift(i))],s);
        roidata2(:,i)=roits{i}([((shift(i)+1):end),(1:shift(i))],t);
    end
    
    nullval=0;
    for i=1:N_subj
        a = roidata1(:,i);
        for j=1:N_subj
            if i~=j
                b = roidata2(:,j);
                nullval = nullval + atanh(sum(sum(a.*b))/div); 
            end
        end
    end
    
    nulldist(iter)=nullval/N_subj_pairs;
    
end

%fprintf('Computing thresholds\n')
%[Th,Th_info,pvals_Th,pvals]=compute_pvals_kauppi(nulldist,mean_conn,0);

end

