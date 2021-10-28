function f = testnetboundary(net,maxnodes,coord)

% Generate correct mesh
[nodes elems msh gm nBNodes] = getMesh(coord);

% Generate input data based on boundary points
nodeIndexIncrement=1;
for nodeIndex=1:nBNodes
    inp{1}(nodeIndexIncrement)  = scaledata(nodes(1,nodeIndex),0,2);
    inp{1}(nodeIndexIncrement+1)= scaledata(nodes(2,nodeIndex),0,2);
    nodeIndexIncrement=nodeIndexIncrement+2;
end

% Since there were variable number of inputs while training, I'm going to pad the exta
% inputs with zeros
newInput(1:maxnodes*2)=0;
for i = 1:1
    actualSize=size(inp{i},2);
    if actualSize < 2*maxnodes
       newInput(1:actualSize)=inp{i};
       newInput(actualSize+1:maxnodes)=0;
    else
       newInput(:)=inp{i}; 
    end
end

size(newInput)


newOut   = net(newInput');
size(newOut)
x        = [ unscaledata(newInput(1:2:end),0,2) newOut(1:2:end)'];
y        = [ unscaledata(newInput(2:2:end),0,2) newOut(2:2:end)'];

%Eliminate all the unwanted points
k=1;
for i =1:size(x,2)
    if inpolygon(x(i),y(i),unscaledata(newInput(1:2:end),0,2), unscaledata(newInput(2:2:end),0,2))
        pointToAdd.x(k)=x(i);
        pointToAdd.y(k)=y(i);
        k=k+1;
    end
end


% 
 tri = delaunay(pointToAdd.x,pointToAdd.y);
 triplot(tri,pointToAdd.x,pointToAdd.y)
 axis equal


% plot(x,y,'.');
% 
% tri = delaunay(x,y);
% triplot(tri,x,y)
% axis equal
