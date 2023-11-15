function [inReg,dist] = distInRegion(source,dest,polygon)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    numPoints = 1000; % You can adjust the number of points for higher accuracy
    xPoints = linspace(source(1), dest(1), numPoints);
    yPoints = linspace(source(2), dest(2), numPoints);
    
    dist=0;
    % Check if any point on the line is inside the polygon
    pointsInside = inpolygon(xPoints, yPoints, polygon(:,1), polygon(:,2));
    
    % Find the indices of points inside the polygon
    indicesInside = find(pointsInside);
    
    % Check if any points are inside the polygon
    if ~isempty(indicesInside)
        inReg=true;
        % Calculate the distance between the first and last points inside the polygon
        distanceFirstLast = sqrt((xPoints(indicesInside(end)) - xPoints(indicesInside(1)))^2 + ...
                                  (yPoints(indicesInside(end)) - yPoints(indicesInside(1)))^2);
    
        % Display the result
        dist=distanceFirstLast;
    else
        inReg=false;
    end
end

