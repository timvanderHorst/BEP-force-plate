function v = vecnorm(vector)
%VECNORM - normalizes each row in a vector
    v = zeros(size(vector,2),1);
    for i = 1 : size(vector,2)
        v(i) = norm(vector(:,i));
    end
end

