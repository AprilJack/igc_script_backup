function rankplot(bcn,prefix,donotplot)
%plot rank plot of counts
%% read in sequences and counts of data split into the different libraries

%bcn=212:229;
for i=1:length(bcn)
data(i).rawcounts=dlmread([prefix,'BC',int2str(bcn(i)),'_counts.txt']);
end
save('rawdata.mat','data','-v7.3');

%use this to look at the sequence rank plot of every libary and choose a threshold for preprocessing.sh 
load rawdata
if ~exist('donotplot')
    figure;for i=1:length(bcn);loglog(data(i).rawcounts);title(int2str(i));pause();end
elseif donotplot~=1
    figure;for i=1:length(bcn);loglog(data(i).rawcounts);title(int2str(i));pause();end
end

