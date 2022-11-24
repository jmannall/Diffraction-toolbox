function NNinput = CreateNNinput(geometry)

    NNinput = ([deg2rad(geometry.wedgeIndex), deg2rad(geometry.bendingAngle), deg2rad(geometry.minAngle), geometry.wedgeLength, geometry.rS, geometry.rR, geometry.zS, geometry.zR])';
end