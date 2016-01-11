clear all 
%Example 1
% %make toy data 
% [startstop,X,y] = LinearToy();
% %Create a KernelRidgeRegressio object 
% ker=['lin'];
% parameters=[0];
% %may have to run multiple times for good results 
% DTAK= DTAKCluster(ker,X,parameters,startstop)
% %Number of clusters
% Clusters=2
% [label, energy] = knkmeans( DTAK,Clusters);
% PlotTime( X ,label,Clusters,startstop)


%Example 2
[startstop,X,y]= NonLinearToy();
ker=['rbf']
parameters=[1000];
%may have to run multiple times for good results 
DTAK= DTAKCluster(ker,X,parameters,startstop)
%Number of clusters
Clusters=2
[label, energy] = knkmeans( DTAK,Clusters);
PlotTime( X ,label,Clusters,startstop)
