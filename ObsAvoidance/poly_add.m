function p = poly_add(p1,p2)
if size(p1,2) == 1
    p1 = p1';
end
if size(p2,2) == 1
    p2 = p2';
end

if length(p1) >= length(p2)
    p = p1 + [zeros(1,length(p1) - length(p2)),p2];
else
    p = [zeros(1,length(p2) - length(p1)),p1] + p2;
end

end