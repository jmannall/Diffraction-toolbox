%% Plot geometry from corner and plane corner data

function PlotGeometry(corners, planeCorners, source, receiver, vSources, intersections, validPath)

    [numPlanes, numCornersPerPlane] = size(planeCorners);
    numCorners = numPlanes * numCornersPerPlane;

    % Remove any trailing zeros from smaller planes
    planeCorners = max(planeCorners, 1);

    % Create planes and normal vectors
    plane = corners(planeCorners',:);
    idx1 = 1:numCornersPerPlane:numCorners;
    idx2 = 2:numCornersPerPlane:numCorners;
    idx3 = 3:numCornersPerPlane:numCorners;
    normals = normr(cross(plane(idx1,:) - plane(idx3,:), plane(idx1,:) - plane(idx2,:)));
    
    % Plot figure
    figure
    plot3(source(:,1), source(:,2), source(:,3), 'o')
    hold on
    plot3(receiver(:,1), receiver(:,2), receiver(:,3), 'o')
    if nargin > 4
        plot3(vSources(:,1), vSources(:,2), vSources(:,3), 'o')
        numVSources = size(vSources, 1);
        for j = 1:numVSources
            if validPath(j)
                intersection = intersections{j};
                vSourcesPlot = [receiver; vSources(j,:)];
                pathPlot = [source; intersection; receiver];
                plot3(vSourcesPlot(:,1), vSourcesPlot(:,2), vSourcesPlot(:,3), 'r')
                plot3(pathPlot(:,1), pathPlot(:,2), pathPlot(:,3), 'g')
            end
        end
    end
    legend('Source', 'Receiver')
    % Loops through planes
    for j = 1:numPlanes
        i = (j - 1) * numCornersPerPlane + 1;
        planePlot = [plane(i:i + numCornersPerPlane - 1,:); plane(i,:)];
        planeCentre = [sum(planePlot(1:numCornersPerPlane,1)) / numCornersPerPlane, sum(planePlot(1:numCornersPerPlane,2)) / numCornersPerPlane, sum(planePlot(1:numCornersPerPlane,3)) / numCornersPerPlane];
        normalsPlot = [planeCentre; planeCentre - normals(j,:)];
        plot3(planePlot(:,1), planePlot(:,2), planePlot(:,3))
        plot3(normalsPlot(:,1), normalsPlot(:,2), normalsPlot(:,3), 'k')
    end
%     xlim([min(plane(:,1)) - 2 max(plane(:,1)) + 2])
%     ylim([min(plane(:,2)) - 2 max(plane(:,2)) + 2])
%     zlim([min(plane(:,3)) - 2 max(plane(:,3)) + 2])
    grid on
    hold off
end