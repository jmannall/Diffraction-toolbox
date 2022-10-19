%% Edge geometry

function [ri, rd, zi, zd, li, ld, phii, apex, inedge, zedgelo, zedgehi, thetai, thetad, v] = EDedgeparams(verts, ps, pk, edgedata, k)
    ev = verts(2,:) - verts(1,:);
    [ui, dzi] = EDupright(verts(1,:), verts(2,:), ps, true);
    [ud, dzd] = EDupright(verts(1,:), verts(2,:), pk, true);
    ri = norm(ui);
    rd = norm(ud);
    dz = abs(dzd - dzi);
    if (ri < 1e-7 || rd < 1e-7)
        msg = 'Source/sink almost on edge';
        error(msg)
    end
    zi = -dz / ((rd / ri) + 1);     % zi < 0
    zd = dz + zi;
    li = sqrt(zi * zi + ri * ri);
    ld = sqrt(zd * zd + rd * rd);
    phii = asin(ri / li);
    if sign(dzd - dzi) < 0
        coordrelsign = -1;
    else
        coordrelsign = 1;
    end
    dz0 = dzi - coordrelsign * zi;
    dz2 = norm(ev);
    apex = EDnormalise(ev) * dz0 + verts(1,:);
    inedge = (0 <= dz0) && (dz0 <= dz2);
    zedge1 = -dz0 * coordrelsign;
    zedge2 = (dz2 - dz0) * coordrelsign;
    zedgelo = min(zedge1, zedge2);
    zedgehi = max(zedge1, zedge2);
    % Angles
    thetaw = 2 * pi  - edgedata.closwedangvec(k);
    v = pi / thetaw;
    check = edgedata.edgenvecs(k,:);
    planeline = EDnormalise(cross(check, ev));
    thetai = EDacos3(planeline, check, EDnormalise(ui));
    thetad = EDacos3(planeline, check, EDnormalise(ud));
    if (thetai > thetad)
        % Angle taken from wrong plane
        thetad = thetaw - thetad;
        thetai = thetaw - thetai;
    end
end