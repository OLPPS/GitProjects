function ran = range(vec)
[Mi, Ma] = minmax(vec);
ran = Ma - Mi;
end

function [Mi, Ma] = minmax(vec)

Mi = min(vec);
Ma = max(vec);

end