%
% Copyright (c) 2016, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YTEA101
% Project Title: Particle Swarm Optimization Video Tutorial
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer and Instructor: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function z = FilterError(x, tfvalue, fs)
    
    [tf, fvec] = IIRFilter([x(1); x(2)], [x(3); x(4)], x(5), fs);
    z = Error(tf, tfvalue, fvec);

end
