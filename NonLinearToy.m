function [startstop,X,y] = NonLinearToy();

%this function creates two classes or clusters of timer series each class %each class or cluster has approximately  1000 samples

% X:Column consist of a different time series dimension the index of the individual time series is determined by    “StopStart”.
%StartStop; (number of time series x  2)  :array, the first column contains the starting address of the time series and the second address contains the stopping address in the time series.
%For example for the series k  StartStop [k,2]=  StartStop [k+1,1]+1 . i.e the  stopping address of series k is one less then  the starting address of series k+1
%y: lables of each time series  

%Empty Feature matrix
X=zeros(2000,2);
%lables of seqence
y=[];
%stop and starting index of each time seresi
STOP=[];
START=[];
%Initialize total size of all sequence
TotalSize=0;
%Starting address of first sequence
START(1)=1;

n=0;
while TotalSize<1000    
    
    %iteration of each sequence 
    n=n+1;
    
    %Number of samples of each time series
    NumberSamples=round(rand(1)*100)
    
    % total length of the n-th sequences
    TotalSize=TotalSize+NumberSamples;
    % stoping address of seqence 
    STOP(n)=START(n)+NumberSamples;
    
    % starting index of next seqence 
    START(n+1)=STOP(n)+1;
    
    %seqwence lable 
    y(n)=0;
    
    %load  array with toy data generated using a linear function  
    
    r = 50; % Radius
    t = 2*pi*[0:1:NumberSamples]/NumberSamples;  % Angle

    
    X(START(n):STOP(n),1)=r.*cos(t)+5*randn(1,1)+100;
    X(START(n):STOP(n),2)=r.*sin(t)+5*randn(1,1)+100;
  
    
    
end



while TotalSize<2000
    %iteration of each sequence 
    n=n+1;
    
    %Number of samples of each time series
    NumberSamples=round(rand(1)*100)
    
    % total length of the n-th sequences
    TotalSize=TotalSize+NumberSamples;
    
    % stoping address of seqence 
    STOP(n)=START(n)+NumberSamples;
    
     % starting index of next seqence 
    START(n+1)=STOP(n)+1;
    
     %seqwence lable 
    y(n)=1;
    
    r = 100; % Radius
    t = 2*pi*[0:1:NumberSamples]/NumberSamples;  % Angle

    
    X(START(n):STOP(n),1)=r.*cos(t)+10*randn(1,1);
    X(START(n):STOP(n),2)=r.*sin(t)+10*randn(1,1);
  
    
    
    
    
end

startstop=zeros(length(STOP(1:end-1)),2);

startstop(:,1)=START(1:end-2);
startstop(:,2)=STOP(1:end-1);




for i = 1:length(STOP)-1    %# Loop 6 times
    
    hold on
    if y(i)==0
    
        plot(X(startstop(i,1):startstop(i,2),2),X(startstop(i,1):startstop(i,2),1),'o');  %# Plot each column with a
    
    else 
        plot(X(startstop(i,1):startstop(i,2),2),X(startstop(i,1):startstop(i,2),1),'or');
    end
  
end

title('Toy Data:True Non-Linear Centroids ' )
