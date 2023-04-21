%% Plot geometry from corner and plane corner data

function PlotIS(room, op)

    vargin = fieldnames(room);
    nargin = length(vargin);

    figure

    %% Source receiver
    
    if op.sr
        s = room.source;
        r = room.receiver;
    
        plot3(s(:,1), s(:,2), s(:,3), 'o')
        hold on
        plot3(r(:,1), r(:,2), r(:,3), 'o')
        legendText = {'Source', 'Receiver'};
    end

    %% Geometry

    if op.geometry
        corners = room.corners;
        planeCorners = room.planeCorners;
    
        [numPlanes, numCornersPerPlane] = size(planeCorners);
        numCorners = numPlanes * numCornersPerPlane;
    
        % Remove any trailing zeros from smaller planes
        %planeCorners = max(planeCorners, 1);
    
        % Create planes and normal vectors
        plane = corners(max(planeCorners, 1)',:);
        idx1 = 1:numCornersPerPlane:numCorners;
        idx2 = 2:numCornersPerPlane:numCorners;
        idx3 = 3:numCornersPerPlane:numCorners;
        normals = normr(cross(plane(idx1,:) - plane(idx3,:), plane(idx1,:) - plane(idx2,:)));
    
        % Loops through planes
        for j = 1:numPlanes
            i = (j - 1) * numCornersPerPlane + 1;
            numCornersInPlane = find(planeCorners(j,:), 1, 'last');
            planePlot = [plane(i:i + numCornersInPlane - 1,:); plane(i,:)];
            planeCentre = mean(planePlot);
            normalsPlot = [planeCentre; planeCentre - normals(j,:)];
            plot3(planePlot(:,1), planePlot(:,2), planePlot(:,3), 'k', 'LineWidth', 1.5)
            plot3(normalsPlot(:,1), normalsPlot(:,2), normalsPlot(:,3), 'k')
            if j == 1
                legendText = [legendText, {'Planes', 'Normals'}];
            else
                legendText = [legendText, {'', ''}];
            end
        end
        xlim([min(plane(:,1)) - 2 max(plane(:,1)) + 2])
        ylim([min(plane(:,2)) - 2 max(plane(:,2)) + 2])
        zlim([min(plane(:,3)) - 2 max(plane(:,3)) + 2])
        grid on
    end

    %% Direct

    if op.direct
        dir = [s; r];
    
        if room.lineOfSight
            plot3(dir(:,1), dir(:,2), dir(:,3), 'b')
        else
            plot3(dir(:,1), dir(:,2), dir(:,3), '--b')
        end
        legendText = [legendText, {'Direct'}];
    end
    
    if op.diff && isfield(room, 'diff')
        numPaths = length(room.diff);
        for i = 1:numPaths
            plot3(room.diff{i}(:,1), room.diff{i}(:,2), room.diff{i}(:,3), 'r')
            if i == 1
                legendText = [legendText, {'1st-order diff'}];
            else
                legendText = [legendText, {''}];
            end
        end
    end

    if op.hodDiff && isfield(room, 'hodDiff')
        diffOrder = length(room.hodDiff);
        for j = 2:diffOrder
            numPaths = length(room.hodDiff{j});
            for i = 1:numPaths
                plot3(room.hodDiff{j}{i}(:,1), room.hodDiff{j}{i}(:,2), room.hodDiff{j}{i}(:,3), 'Color', [0.6 0 0])
                if i == 1
                    switch j
                        case 2
                            text = 'nd';
                        case 3
                            text = 'rd';
                        otherwise
                            text = 'th';
                    end
                    legendText = [legendText, {[num2str(j), text, '-order diff']}];
                else
                    legendText = [legendText, {''}];
                end
            end
        end
    end



    legend(legendText);

    


    if op.reflect
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
    direct = [source; receiver];
    % Plot figure
    figure
    plot3(source(:,1), source(:,2), source(:,3), 'o')
    hold on
    plot3(receiver(:,1), receiver(:,2), receiver(:,3), 'o')
    legend('Source', 'Receiver')
    if nargin > 4 && lineOfSight
        plot3(direct(:,1), direct(:,2), direct(:,3), 'b')
    end
    if nargin > 5
        plot3(vSources(:,1), vSources(:,2), vSources(:,3), 'o', 'Color', 'm')
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
    if nargin > 8
        numEdgePaths = length(edgePaths);
        for j = 1:numEdgePaths
            plot3(edgePaths{j}(:,1), edgePaths{j}(:,2), edgePaths{j}(:,3), 'b')
        end
    end
    % Loops through planes
    for j = 1:numPlanes
        i = (j - 1) * numCornersPerPlane + 1;
        planePlot = [plane(i:i + numCornersPerPlane - 1,:); plane(i,:)];
        planeCentre = [sum(planePlot(1:numCornersPerPlane,1)) / numCornersPerPlane, sum(planePlot(1:numCornersPerPlane,2)) / numCornersPerPlane, sum(planePlot(1:numCornersPerPlane,3)) / numCornersPerPlane];
        normalsPlot = [planeCentre; planeCentre - normals(j,:)];
        plot3(planePlot(:,1), planePlot(:,2), planePlot(:,3))
        plot3(normalsPlot(:,1), normalsPlot(:,2), normalsPlot(:,3), 'k')
    end
    xlim([min(plane(:,1)) - 2 max(plane(:,1)) + 2])
    ylim([min(plane(:,2)) - 2 max(plane(:,2)) + 2])
    zlim([min(plane(:,3)) - 2 max(plane(:,3)) + 2])
    grid on
    hold off
    end
end