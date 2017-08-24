vec2mat(v::AbstractVector) = reshape(v, (length(v),1))
vec2mat(T::Type, v::AbstractVector) = vec2mat(T.(v))
