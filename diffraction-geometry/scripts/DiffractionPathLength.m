function pathLength = DiffractionPathLength(receivers, source, apex)

    pathLength = PathLength(receivers, apex) + PathLength(apex, source);
end