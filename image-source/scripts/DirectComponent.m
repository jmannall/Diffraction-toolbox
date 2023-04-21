function [direct, plot] = DirectComponent(room, plot)

    i = 0;
    direct.valid = ~CheckForObstruction(room.source, room.receiver, room, i);
    plot.lineOfSight = direct.valid;
    
    if direct.valid
        direct.pathLength = PathLength(room.receiver, room.source);
    end
end