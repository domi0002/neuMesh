clear all
clc
k=1;
maxsize=0;
% Generate input and output dataset

dmin=0.5;
dmax=2.0;
maxBoundaryNodes=0;

for i = dmin:0.5:dmax
%     inp(k,1)=0;
%     inp(k,2)=1;
%     inp(k,3)=i;
%     inp(k,4)=0;
%     inp(k,5)=0;
%     inp(k,6)=0;
%     inp(k,7)=1;
%     inp(k,8)=1;
%     singleinput(k)=i;
     % scale the input 
%     singleinput(k)=-1+2*(i-dmin)/(dmax-dmin);
  
    [nodes elems msh gm nBNodes] = getMesh(i);
     nodeIndexIncrement=1;
     for nodeIndex=1:nBNodes
%              inp(k,nodeIndexIncrement)=nodes(1,nodeIndex);
%              inp(k,nodeIndexIncrement+1)=nodes(2,nodeIndex);
              inp{k}(nodeIndexIncrement)= scaledata(nodes(1,nodeIndex),0,dmax);
              inp{k}(nodeIndexIncrement+1)=scaledata(nodes(2,nodeIndex),0,dmax);
       
             nodeIndexIncrement=nodeIndexIncrement+2;
     end
     
     %Determine max number of boundary nodes in all sets
     maxBoundaryNodes=max(maxBoundaryNodes,nBNodes);
    
    % Truncated outputs, i.e. points excluding Bpoints 
    zi=1;
    for z = nBNodes+1:size(nodes,2)
        
        %zi=2*z-1;
        
        out{k}(zi)   = nodes(1,z);
        out{k}(zi+1) = nodes(2,z);
        %el{k}        = elems;
        
        zi=zi+2;
        
    end
    
    % These are the full outputs- will be used later on if required
    zi=1;
    for z = 1:size(nodes,2)
        
        %zi=2*z-1;
        fullOut{k}(zi)   = nodes(1,z);
        fullOut{k}(zi+1) = nodes(2,z);
        el{k}        = elems;
        zi=zi+2;
    end
    
    
    
    %out{k} = [nodes(1,:) nodes(2,:)];
    %out(k,:)=out1;
    
    if size(out{k},2) > maxsize
        maxsize = size(out{k},2);
    end
        
    k=k+1;
end
kmax=k-1;

% Since there are variable number of inputs, I'm going to pad the exta
% inputs with zeros
newInput(1:kmax,1:maxBoundaryNodes*2)=0;
for i = 1:kmax  
    actualSize=size(inp{i},2);
    if actualSize < 2*maxBoundaryNodes
       newInput(i,1:actualSize)=inp{i};
       newInput(i,actualSize+1:maxBoundaryNodes)=0;
    else
       newInput(i,:)=inp{i}; 
    end
end



newOutput(1:kmax,1:maxsize)=0;
% Pad with zeros?
for i = 1:kmax
    actualSize= size(out{i},2);
    if actualSize < maxsize
        nElementsToPad=maxsize-actualSize;
        newOutput(i,1:actualSize)=out{i};
         for k = 1:nElementsToPad/2
          %for k = size(el{i},2):-1:size(el{i},2)-nElementsToPad/2
            % Pick these vertices vTp
            vTp = el{i}(:,k);
            newOutput(i,actualSize+2*k-1)= (1/3)*( fullOut{i}(2*vTp(1)-1) +... 
                                                   fullOut{i}(2*vTp(2)-1) +...
                                                   fullOut{i}(2*vTp(3)-1));
            newOutput(i,actualSize+2*k)  = (1/3)*( fullOut{i}(2*vTp(1)) +... 
                                                   fullOut{i}(2*vTp(2)) +...
                                                   fullOut{i}(2*vTp(3)));
         end
%        newOutput(i,actualSize+1:end)=10*ones(1,nElementsToPad);
    else
        newOutput(i,:)=out{i};
    end
end

%  x = newOutput(3,1:2:end);
%  y = newOutput(3,2:2:end);
%   tri = delaunay(x,y);
%  triplot(tri,x,y);
%  axis equal


% Train data here (lower neurons seems better!)
net=feedforwardnet(20);
net.trainFcn='trainlm';
net.trainParam.goal=1e-9;
net.trainParam.epochs=1000;
net.trainParam.min_grad=1e-9;
net.trainParam.max_fail=50;
net.divideParam.trainRatio=1;
net.divideParam.valRatio=0;
net.divideParam.testRatio=0;
%net.performFcn='msereg'; %net.performParam.ratio=0.5;
%[net, tr]= train(net, singleinput, newOutput');
[net, tr]= train(net, newInput', newOutput');

%desiredInput=0.51;
%scaledDesiredInput=-1+2*(desiredInput-dmin)/(dmax-dmin)
%testnet(net,scaledDesiredInput)

desiredInput=1.35;
testnetboundary(net,maxBoundaryNodes,desiredInput);

    
    