function gr = gram_row(i, s, n, quad_form)
gr = zeros([n 1], class(s));
point1 = s(:,i);
for j = 1:i
    point2 = s(:,j);
    gr(i) = -sum([point1 point2]'.*(quad_form*[point1 point2]), 'all');
end
end