function fn = getMesh(coordinate)
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
pdegplot(gm, 'EdgeLabels','on')
axis equal

generateMesh(pdem);

pdeplot(pdem)
