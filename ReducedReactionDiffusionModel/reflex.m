function [i1, j1] = reflex(a, b, R, H, h)
syms r
curve = @(r) -(H/R)*r+H;
d = diff(curve, r);
eqn = a-r + b*d -curve*d == 0;
S = vpasolve(eqn,r,[0 R]);
r0 = double(S);
z0 = curve(r0);
a1 = 2*r0-a;
b1 = 2*z0-b;
i1 = round(a1/h);
j1 = round(b1/h);
% disp([i1,j1])
if curve(i1*h)<=j1*h
    if curve(i1*h)<=(j1-1)*h
        i1 = i1-1;
    elseif curve((i1-1)*h)<=j1*h
        j1 = j1-1;
    else
        if a1-(i1-1)*h>=b1-(j1-1)*h
            i1 = i1-1;
        else
            j1 = j1-1;
        end
    end
end
if i1 == R/h
    disp('y');
    i1 = i1-1;
end
if j1<=0
    j1=0;
end
