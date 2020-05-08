# quantum time evolution

using FFTW

mutable struct Operators
    V::Vector{Float64}
    R::Vector{Complex{Float64}}
    K::Vector{Complex{Float64}}
    Operators(L) = new(zeros(L), zeros(L), zeros(L))
end

function oprinit(V::Function)
    opr::Operators = Operators(res)
    opr.V::Vector{Float64} .+= V.(xspc)
    opr.R::Vector{Complex{Float64}} .+= exp.(-im .* opr.V .* 0.5(imtm ? dt/im : dt))
    opr.K::Vector{Complex{Float64}} .+= exp.(-im .* 0.5(kspc.^2) .* (imtm ? dt/im : dt))
    opr
end

function energy(wfn::Vector{Complex{Float64}}, opr::Operators)
    wfn_k::Vector{Complex{Float64}} = fft(wfn)
    wfn_c::Vector{Complex{Float64}} = conj.(wfn)
    Ek::Vector{Complex{Float64}} = wfn_c .* ifft(0.5(kspc.^2) .* wfn_k)
    Er::Vector{Complex{Float64}} = wfn_c .* opr.V .* wfn
    E::Float64 = 0.
    for i = 1:res
        E += hypot(Ek[i] + Er[i])
    end
    E / res
end

function anim(wfn::Vector{Complex{Float64}}, opr::Operators)
    for i = 1:(res - 1)
        println("l3 $((i-mid)*dx) $(real(wfn[i])) $(imag(wfn[i])) $(((i+1)-mid)*dx) $(real(wfn[i+1])) $(imag(wfn[i+1]))")
        println("l3 $((i-mid)*dx) $(opr.V[i]) 0.0 $(((i+1)-mid)*dx) $(opr.V[i+1]) 0.0")
    end
    if enrg
        println("T -0.8 0.8")
        println("E = $(energy(wfn, opr))")
    end
    println("F")
end

function split!(wfn::Vector{Complex{Float64}}, opr::Operators)
    wfn::Vector{Complex{Float64}} .*= opr.R
    wfnk::Vector{Complex{Float64}} = fft(wfn)
    wfnk::Vector{Complex{Float64}} .*= opr.K
    wfn::Vector{Complex{Float64}} .= ifft(wfnk)
    wfn::Vector{Complex{Float64}} .*= opr.R
end

function spawn!(wfn::Vector{Complex{Float64}})
    I = div(res, 3)
    for i = 0:λ
        wfn[i + I] = sinpkt(i)
        # wfn[i + I] = exppkt(i)
            # sin(π * i / λ) * exp(im * p^2 * π * i / 4 / λ)
    end
end

function normalize!(wfn::Vector{Complex{Float64}})
    sum::Float64 = 0.
    for i = 1:res
        sum += hypot(wfn[i])
    end
    wfn ./= sum * dx
    wfn
end

sinpkt(x) = sin(π * x / λ) * exp(-im * p^2 * π * x / λ)
exppkt(x) = A * exp(-((x - λ/2)/σ)^2) / (σ * sqrt(2π))


const res  = 2^8::Int64
const len  = 20.::Float64
const mid  = div(res, 2)
const dx   = len / res
const dk   = 2π / len
const dt   = 1e-3::Float64
const xspc = Vector{Float64}((-mid : mid - 1) * dx)
const kspc = Vector{Float64}(vcat(0 : mid - 1, -mid:-1) * dk)

const k = 1.5
const p = 3.0::Float64
const λ = div(res, 6)
const σ = λ/6
const A = 20

const enrg = true
const imtm = false

function main()
    V(x) = 20x^2
    # V(x) = begin x <= (mid + 2)
    wfn::Vector{Complex{Float64}} = zeros(res)
    spawn!(wfn)
    while true
        anim(wfn, opr)
        split!(wfn, opr)
    end
end

@time main()
