function [fn1, fn2, fn3, fn4, fn5 ] = getMesh(coordinate)
% Create Training data set

R1 = [3,4,0,1,coordinate,0,0,0,1,1]';
       
geom = R1;
% Names for the two geometric objects
ns = (char('R1'))';
   
% Set formula
sf = 'R1';
% Create geometry
gd = decsg(geom,sf,ns);
% Now create the class geometry representation
pdem = createpde();
gm = geometryFromEdges(pdem, gd)
% Plot the geometry
%pdegplot(gm, 'EdgeLabels','on')
%axis equal

hmax=0.2;
msh=generateMesh(pdem,'Hmax',hmax, 'GeometricOrder','linear');

fn1 =msh.Nodes;
fn2 =msh.Elements;
fn3 = msh;
fn4 = gm;

%pdeplot(pdem)

% This will be the number of boundary nodes
approxBoundaryEdges=(1/hmax) + sqrt((coordinate-1)^2 + 1^2)/hmax + coordinate/hmax  + 1/hmax;
fn5=ceil(approxBoundaryEdges);



