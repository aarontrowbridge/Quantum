# path integral double slit

module Paths

using LinearAlgebra

export propagate

S(xi, xf) = 0.5 * norm(xf - xi)^2

K(xi, xf) = exp(im * S(xi, xf))

source(ds) = (div(ds[1] + 1, 2), 1)

const dims = (11, 20)

function propagate(ds::Tuple{Int64, Int64})

    space = zeros(Complex, ds)

    (xsrc, ysrc) = source(ds)

    space[xsrc, ysrc] = 1.

    for i = 1:size(space)[1]
        xi = [xsrc, ysrc]
        xf = [i - xsrc, 1]
        space[i,1] = K(xi, xf)
    end

    for j = 2:size(space)[2]
        for k = 1:size(space)[1]
            for i = 1:size(space)[1]
                xi = [i - xsrc, j - 1]
                xf = [k - xsrc, j]
                space[k, j] += space[i, j - 1] * K(xi, xf)
            end
        end
        space[:,j-1] ./= sum(abs.(space[:,j-1]))
    end

    space[:,end] ./= sum(abs.(space[:,end]))

    space
end

end
