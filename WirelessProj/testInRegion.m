% Define the vertices of the polygon
polygonVertices = [1, 1; 1, 5; 5, 5; 5, 1];
polygonVertices1 = [5, 5; 5,6; 6,5; 6, 6];

% Define the two points of the line segment
point1 = [7,8];
point2 = [2,2];
disp("start")
% Check if any point on the line segment is inside the polygon
distance=pdist2(point1,point2);


[in1, dist1]=distInRegion(point1,point2,polygonVertices);

[in2, dist2]=distInRegion(point1,point2,polygonVertices1);

disp(distance)
disp(distance-dist1-dist2)

