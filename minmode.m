function [c] = minmode(a,b)
if a*b <= 0
    c = 0;
else 
    c = sign(a)*min(abs(a),abs(b));
end
end