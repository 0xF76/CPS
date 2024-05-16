module CPS

using LinearAlgebra

author = Dict{Symbol, String}(
    :index => "414702",
    :name  => "Jakub Bubak",
    :email => "jakubbubak@student.agh.edu.pl",
    :group => "1",
)

# Sygnały ciągłe
cw_rectangular(t::Real; T=1.0)::Real = abs(t) < T/2 ? 1 : 0
cw_triangle(t::Real; T=1.0)::Real = abs(t) < T ? 1 - 1/T*abs(t) : 0
cw_literka_M(t::Real; T=1.0)::Real = abs(t) < T/2 ? abs(t) : 0
cw_literka_U(t::Real; T=1.0)::Real = missing

ramp_wave(t::Real)::Real = t - floor(t)
sawtooth_wave(t::Real)::Real = -(t-floor(t))
triangular_wave(t::Real)::Real = 4*abs(t-floor(t+1/2))
square_wave(t::Real)::Real =  sign(sin(2*π*t))
pulse_wave(t::Real, ρ::Real=0.2)::Real = (t - floor(t)) ≤ ρ ? 1 : 0
impulse_repeater(g::Function, t1::Real, t2::Real)::Function = x -> g(mod(x - t1, t2 - t1) + t1)

function ramp_wave_bl(t; A=1.0, T=1.0, band=20.0)
    sum = 0
    k = 1
    while (arg = (2π*k)/T) ≤ band * 2π
        sum += (-1)^k * sin(arg * t) / k
        k += 1
    end
    return -2 * A / π * um
end

function sawtooth_wave_bl(t; A=1.0, T=1.0, band=20.0)
    sum = 0
    k = 1
    while (arg = (2π*k)/T) ≤ band * 2π
        sum += (-1)^k * sin(arg * t) / k
        k += 1
    end
    return 2 * A / π * sum
end

function triangular_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

function square_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

function pulse_wave_bl(t; ρ=0.2, A=1.0, T=1.0, band=20.0)
    missing
end


function impuse_repeater_bl(g::Function, t0::Real, t1::Real, band::Real)::Function
    missing
end

function rand_signal_bl(f1::Real, f2::Real)::Function
    missing
end




# Sygnały dyskretne
kronecker(n::Integer)::Real = n == 0 ? 1 : 0
heaviside(n::Integer)::Real = n ≥ 0 ? 1 : 0

# Okna
rect(N::Integer)::AbstractVector{<:Real} = ones(N)
triang(N::Integer)::AbstractVector{<:Real} = [1 - (2abs(n - ((N - 1) / 2))) / (N - 1) for n = 0:N-1]
hanning(N::Integer)::AbstractVector{<:Real} = [0.5*(1 - cos((2π*n)/(N-1))) for n = 0:N-1]
hamming(N::Integer)::AbstractVector{<:Real} = [0.54 - 0.46*cos((2π*n)/(N-1)) for n = 0:N-1]
blackman(N::Integer)::AbstractVector{<:Real} = [0.42 - 0.5*cos((2π*n)/(N-1)) + 0.08*cos((4π*n)/(N-1)) for n = 0:N-1]

# Parametry sygnałów
mean(x::AbstractVector)::Number = sum(x)/length(x)
peak2peak(x::AbstractVector)::Real = maximum(x) - minimum(x)
energy(x::AbstractVector)::Real = sum(x.^2)
power(x::AbstractVector)::Real = energy(x)/length(x)
rms(x::AbstractVector)::Real = √power(x)

function running_mean(x::AbstractVector, M::Integer)::Vector
    missing
end

function running_energy(x::AbstractVector, M::Integer)::Vector
    missing
end

function running_power(x::AbstractVector, M::Integer)::Vector
    missing
end



# Próbkowanie
function interpolate(
    m::AbstractVector,
    s::AbstractVector,
    kernel::Function=sinc
)::Real
    missing
end

# Kwantyzacja
quantize(x::Real; L::AbstractVector)::Real = missing
SQNR(N::Integer)::Real = 1.76 + 6.02*N # 6.02*N [dB] is also correct
SNR(Psignal, Pnoise)::Real = 10*log10(Psignal/Pnoise)


# Obliczanie DFT
function dtft(f::Real; signal::AbstractVector, fs::Real)
    val::ComplexF64 = 0.0
    for n in eachindex(signal)
        dtft_val += signal[n] * exp(-2 * f/fs * n)
    end
    
    return val
end

function dft(x::AbstractVector)::Vector
    N = length(x)
    n = 0:N-1
    k = 0:N-1
    return exp.(-1im * 2 * π/N * k * n') * x
end

function idft(X::AbstractVector)::Vector
   N = length(X)
    n = 0:N-1
    k = 0:N-1
    return exp.(-1im * 2 * π/N * k * n')^(-1) * X
end

function rdft(x::AbstractVector)::Vector
   missing
end

function irdft(X::AbstractVector, N::Integer)::Vector
   missing
end

function fft_radix2_dit_r(x::AbstractVector)::Vector
   missing
end

function ifft_radix2_dif_r(X::AbstractVector)::Vector
   missing
end

# TODO: THIS IS ABSOLUTELY NOT THE BEST WAY TO DO IT, BUT IT WORKS
function fft(x::AbstractVector)::Vector
    if length(x) ≤ 1
        return x
    else
        N = length(x)
        Xᵒ = fft(x[2:2:end])
        Xᵉ = fft(x[1:2:end])
        factors = exp.(-2im * π * (0:(N/2 - 1)) / N)
        return [Xᵉ .+ factors .* Xᵒ; Xᵉ .- factors .* Xᵒ]
    end
end

function ifft(X::AbstractVector)::Vector
    idft(X) # Może da rade lepiej?
end

fftfreq(N::Integer, fs::Real)::Vector = missing
rfftfreq(N::Integer, fs::Real)::Vector = missing
amplitude_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = missing
power_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = missing
psd(x::AbstractVector, w::AbstractVector=rect(length(x)), fs::Real=1.0)::Vector = missing

function periodogram(x::AbstractVector, w::AbstractVector=rect(length(x)), fs::Real=1.0)::Vector
    missing
end



function stft(x::AbstractVector, w::AbstractVector, L::Integer)::Matrix
    missing
end


function istft(X::AbstractMatrix{<:Complex}, w::AbstractVector{<:Real}, L::Integer)::AbstractVector{<:Real}
    missing
end

function conv(f::Vector, g::Vector)::Vector
    y = zeros(length(f)+length(g)-1)
    for m = 1:length(f)
        for n = 1:length(g)
            y[m+n-1] += f[m]*g[n]
        end
    end
end

function fast_conv(f::Vector, g::Vector)::Vector
    N = length(f)
    M = length(g)
    K = N + M - 1
    while(length(g) < K)
        append!(g,0)
    end

    while(length(f) < K)
        append!(f,0)
    end

    g = g[end:-1:1]
    g = [g[end]; g[1:end-1]]

    c = zeros(K)

    for i = 1:K
        c[i] = sum(f.*g)    
        g = [g[end]; g[1:end-1]]
    end

    return c
end

function overlap_add(x::Vector, h::Vector, L::Integer)::Vector
    missing
end

function overlap_save(x::Vector, h::Vector, L::Integer)::Vector
    missing
end

function lti_filter(b::Vector, a::Vector, x::Vector)::Vector
    missing
end

end

