%% UTD model

function ret = kouyoumjianlenatten(dinc, ddiff, k)
    num = exp(-1i  * k * (dinc + ddiff));
    denom = sqrt(dinc * ddiff * (dinc + ddiff));
    ret = num / denom;
end