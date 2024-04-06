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
cw_triangle(t::Real; T=1.0)::Real = abs(t) < T/2 ? 1 - 1/(T/2)*abs(t) : 0
cw_literka_M(t::Real; T=1.0)::Real = missing
cw_literka_U(t::Real; T=1.0)::Real = missing

ramp_wave(t::Real)::Real = missing
sawtooth_wave(t::Real)::Real = -(t-floor(t))
triangular_wave(t::Real)::Real = 2*abs(t-floor(t+1/2))
square_wave(t::Real)::Real =  sign(sin(2*π*t))
pulse_wave(t::Real, ρ::Real=0.2)::Real = (t - floor(t)) ≤ ρ ? 1 : 0
impulse_repeater(g::Function, t1::Real, t2::Real)::Function = missing

function ramp_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

function sawtooth_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
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
rect(N::Integer)::AbstractVector{<:Real} = missing
triang(N::Integer)::AbstractVector{<:Real} = missing
hanning(N::Integer)::AbstractVector{<:Real} = missing
hamming(N::Integer)::AbstractVector{<:Real} = missing
blackman(N::Integer)::AbstractVector{<:Real} = missing

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
SQNR(N::Integer)::Real = missing
SNR(Psignal, Pnoise)::Real = Psignal/Pnoise


# Obliczanie DFT
function dtft(f::Real; signal::AbstractVector, fs::Real)
   missing
end

function dft(x::AbstractVector)::Vector
    missing
end

function idft(X::AbstractVector)::Vector
   missing
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

function fft(x::AbstractVector)::Vector
    dft(x) # Może da rade lepiej?
end

function ifft(X::AbstractVector)::Vector
    idft(X) # Może da rade lepiej?
end

end

