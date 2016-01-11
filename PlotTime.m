function [] = PlotTime( X ,label,Clusters,startstop)
%Plots  different time series with labels according to cluster labels 
%Input
%X: data matrix with training samples in rows and
%features in columns, to get the idex of a  particular time series see startstop
%startstop: starting a stoping point of a time seres
%label: labels of each time series 
%Clusters

y=X'


figure
cmap = hsv(Clusters)

for k=1:Clusters

dummy=startstop(label==k,:);
%number of samples that belong to cluster k 

Nsamples=sum(label==k)
hold on
for n=1:Nsamples

plot(y(2,dummy(n,1):dummy(n,2)),y(1,dummy(n,1):dummy(n,2)),'o','Color',cmap(k,:),'MarkerSize',8); 
end
end
title('Cluster Lables') 
hold off

end