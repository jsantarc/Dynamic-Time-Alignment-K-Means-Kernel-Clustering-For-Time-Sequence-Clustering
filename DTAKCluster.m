classdef DTAKCluster
    
    %This objected preforms dynamic time aliment kernel clustering. Declare the object by inputting the set of timer series, timer series index,
    %the kernel you want and the parameters.Then input the kernal matrix K and the number of clusters you want into the method knkmeans(K,numberofclusters)
    %and the data set will be labled
   % by Joseph Santarcangelo 
    
    %Thanks%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % I used code or bits of code form the following people, who I would like to thank:
    %Michael Chen: clustering algothim (sth4nth@gmail.com).
    
    %Quan Wang:Used parts of his dynamic time warping algorithm    
    %(www.mathworks.com/matlabcentral/profile/authors/2797017-quan-wang)
    
    % Gustavo Camps-Valls and Jordi (jordi@uv.es) for the fast kernal calculation 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    properties
        
        %kennel matrix
        K;
        
    end
    
    methods
        function DTAKCluster= DTAKCluster(ker,X,parameters,startstop)
            %ker: 'lin','poly','rbf'
            %X: data matrix with training samples in rows and
            %features in columns, to get the idex of a  particular
            % time series see startstop
            %startstop: starting a stoping point of a time seres
            
            %parameters:
            %If RBF kernel or sam kernel must be one value
            %sigma: width of the RBF kernel
            
            % for polynomial kernel 'poly' the paramters are in the
            % form  parameters=[b a]
            %       b:     bias in the linear and polinomial kernel
            %       d:     degree in the polynomial kernel
            
            %Linear kernels 'lin' parameters=b;
            %       b:     bias in the linear and polinomial kernel
            
            
            %Target:row vector of continuous values to be predicted
            
            % %K:Gram matrix, each element corresponds to the kernel of two feature vectors
            
            % Check stuff
            if (strcmp(ker,'rbf')||strcmp(ker,'sam'))
                
                if  length(parameters)~=1
                    error('Error: RBF kernel and Sam kennel only needs one parameter')
                end
                
                b=1;
                d=1;
                sigma=parameters;
                DTWKernelMatrix=KernelMatrix(X,startstop,b ,d,sigma,ker );
            end
            
            if strcmp(ker,'poly')
                if  length(parameters)~=2
                    error('Error: polynomial kernels need two parameters')
                end
                
                sigma=1;
                b=parameters(1);
                d=parameters(2);
                DTWKernelMatrix=KernelMatrix(X,startstop,b ,d,sigma,ker );
            end
            
            if strcmp(ker,'lin')
                if  length(parameters)~=1
                    error('Error: Linear kernels only need one parameter')
                end
                
                sigma=1;
                d=1;
                b=parameters(1);
                DTWKernelMatrix=KernelMatrix(X,startstop,b ,d,sigma,ker );
            end
            
            DTAKCluster.K=DTWKernelMatrix;
        end
        
        function [label, energy] = knkmeans(DTAKCluster,init)
            % Perform kernel k-means clustering.
            % DTAKCluster:DTAKCluster clustering object  
            % init: k (1 x 1) or label (1 x n, 1<=label(i)<=k)
            % Reference: [1] Kernel Methods for Pattern Analysis
            % by John Shawe-Taylor, Nello Cristianini
            % Written by Michael Chen (sth4nth@gmail.com).
           
            %kernel matrix
            K=DTAKCluster.K;
            
            n = size(K,1);
            if length(init) == 1
                label = ceil(init*rand(1,n));
            elseif size(init,1) == 1 && size(init,2) == n
                label = init;
            else
                error('ERROR: init is not valid.');
            end
            last = 0;
            while any(label ~= last)
                [u,~,label] = unique(label);   % remove empty clusters
                k = length(u);
                E = sparse(label,1:n,1,k,n,n);
                E = bsxfun(@times,E,1./sum(E,2));
                T = E*K;
                Z = repmat(diag(T*E'),1,n)-2*T;
                last = label;
                [val, label] = min(Z,[],1);
            end
            [~,~,label] = unique(label);   % remove empty clusters
            energy = sum(val)+trace(K);
        end
    end
    
    
end


function Km=KernelMatrix(y,startstop,Bias ,Power,sigma,ker )

%Output
%Km:Kernel matrix
%Input
%y:Column consist of a different time series dimension the index of the individual time series is determined by  “StopStart”.

%StartStop; (number of time series x  2)  :array, the first column contains the starting address of the time series and the
%second address contains the stopping address in the time series.
%For example for the series k  StartStop [k,2]=  StartStop [k+1,1]+1 .
%i.e the  stopping address of series k is one less then  the starting address of series k+1
%ker: name of Kernel
%Kernel parameters
%Bias ,Power,sigma,

y=y';


NumberSequences=length(startstop(:,1)) ;

Km=zeros(NumberSequences,NumberSequences);

for i=1:NumberSequences
    
    for j=i:NumberSequences
        
        %find Dynamic Time-Alignment Kernel
        [Kernel,Dist,D,k,w,d]=kerneldtw(y(:,startstop(i,1):startstop(i,2)),y(:,startstop(j,1):startstop(j,2)),Bias ,Power,sigma,ker );
        
        Km(i,j)=Kernel;
        
        Km(j,i)=Km(i,j);
        
    end
end

end

function [Kernel,Dist,D,k,w,d]=kerneldtw(X,X2,Bias ,Power,sigma,ker )
% This fuction caculates the Dynamic Time-Alignment Kernel between two
% sereises modified from Quan Wang
%Input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%X: series one order is mixed up i.e columns time index rows are dimension
%X: series two order is mixed up i.e columns time index rows are dimension
%ker: Kernel type
%Different kennel parameters
%Bias
%Power
%sigma

%Output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Dist is unnormalized distance between t and r
%D is the accumulated distance matrix
%k is the normalizing factor
%w is the optimal path
%d:Inner product matrix
%t is the vector you are testing against
%r is the vector you are testing
%Power:degree in the polynomial kernel
%Bias:bias in the linear and polinomial kernel
%Modification so two detentions can be used
%Get first detentions


[rows,N]=size(X);
[rows,M]=size(X2);
if M==0;
    
    M=1;
end

d=zeros(N,M);
switch ker
    case 'lin'
        d=(X' * X2 + Bias);
        
        %Polynomial kernel
    case 'poly'
        d=(X' * X2 + Bias).^Power;
        
    case 'rbf'
        
        n1sq = sum(X.^2,1);
        n1 = size(X,2);
        n2sq = sum(X2.^2,1);
        n2 = size(X2,2);
        D1 = (ones(n2,1)*n1sq)' + ones(n1,1)*n2sq -2*X'*X2;
        d = exp(-D1/(2*sigma^2));
        
end

%d=(repmat(t(:),1,M)-repmat(r(:)',N,1)).^2; %this replaces the nested for loops from above Thanks Georg Schmitz
D=zeros(size(d));
D(1,1)=d(1,1);

for n=2:N
    D(n,1)=d(n,1)+D(n-1,1);
end
for m=2:M
    D(1,m)=d(1,m)+D(1,m-1);
end
for n=2:N
    for m=2:M
        D(n,m)=d(n,m)+max([D(n-1,m),D(n-1,m-1),D(n,m-1)]);
    end
end

Dist=D(N,M);
n=N;
m=M;
k=1;
w=[];
w(1,:)=[N,M];
while ((n+m)~=2)
    if (n-1)==0
        m=m-1;
    elseif (m-1)==0
        n=n-1;
    else
        [values,number]=max([D(n-1,m),D(n,m-1),D(n-1,m-1)]);
        switch number
            case 1
                n=n-1;
            case 2
                m=m-1;
            case 3
                n=n-1;
                m=m-1;
        end
    end
    k=k+1;
    w=cat(1,w,[n,m]);
end

%Path of max correlation
wnew=fliplr(w')';
Kernel=0;
%Path Normalizing Factor
NormalizingF=0;
L=length(wnew(:,1));

%find the max path
for l=1:L
    
    Kernel=Kernel+d(wnew(l,1),wnew(l,2))*D(wnew(l,1),wnew(l,2));
    NormalizingF=NormalizingF+D(wnew(l,1),wnew(l,2));
end

if NormalizingF==0
    NormalizingF=1;
end

%Normalize kernel
Kernel=Kernel/NormalizingF;
end

function [label, energy] = knkmeans(K,init)
% Perform kernel k-means clustering.
%   K: kernel matrix
%   init: k (1 x 1) or label (1 x n, 1<=label(i)<=k)
% Reference: [1] Kernel Methods for Pattern Analysis
% by John Shawe-Taylor, Nello Cristianini
% Written by Michael Chen (sth4nth@gmail.com).
n = size(K,1);
if length(init) == 1
    label = ceil(init*rand(1,n));
elseif size(init,1) == 1 && size(init,2) == n
    label = init;
else
    error('ERROR: init is not valid.');
end
last = 0;
while any(label ~= last)
    [u,~,label] = unique(label);   % remove empty clusters
    k = length(u);
    E = sparse(label,1:n,1,k,n,n);
    E = bsxfun(@times,E,1./sum(E,2));
    T = E*K;
    Z = repmat(diag(T*E'),1,n)-2*T;
    last = label;
    [val, label] = min(Z,[],1);
end
[~,~,label] = unique(label);   % remove empty clusters
energy = sum(val)+trace(K);
end