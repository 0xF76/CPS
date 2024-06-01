module CPS

using LinearAlgebra
using OffsetArrays
using QuadGK

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
cw_literka_U(t::Real; T=1.0)::Real = abs(t) < T/2 ? t^2 : 0

ramp_wave(t::Real)::Real = t - floor(t)
sawtooth_wave(t::Real)::Real = -2*(t-floor(t+1/2))
triangular_wave(t::Real)::Real = 4*abs((t+1/4)-floor((t+1/4)+1/2)) - 1
square_wave(t::Real)::Real =  sign(sin(2*π*t))
pulse_wave(t::Real, ρ::Real=0.2)::Real = (t - floor(t)) ≤ ρ ? 1 : 0
impulse_repeater(g::Function, t1::Real, t2::Real)::Function = x -> g(mod(x - t1, t2 - t1))

function ramp_wave_bl(t; A=1.0, T=1.0, band=20.0)
    sum = 0
    k = 1
    while (arg = (2π*k)/T) ≤ band * 2π
        sum += (-1)^k * sin(arg * t) / k
        k += 1
    end
    return -2 * A / π * sum
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
    sum = 0
    k = 1
    while (arg = (2π*(2*k-1))) ≤ band * 2π
        sum += ((-1)^k)/(2k-1)^2 * sin(arg*t)
        k+=1
    end
    
    return -8 / π^2 * sum
end

function square_wave_bl(t; A=1.0, T=1.0, band=20.0)
    signal = 0
    n = 1
    while (arg = 2π * (2n - 1) * (1 / T)) < band * 2π
        signal += sin.(arg * t) / (2n - 1)
        n += 1
    end
    signal *= 4 * A / π
    return signal
end

function pulse_wave_bl(t; ρ=0.2, A=1.0, T=1.0, band=20.0)
    sum = 0
    k = 1
    while (arg = 2*π*k/T) ≤ band * 2π
        sum += 1/k * sin(π*k*ρ/T)*cos(arg*t)
        k+=1
    end

    return A*ρ/T + 2*A/π * sum

end


function impulse_repeater_bl(g::Function, t1::Real, t2::Real, band::Real)::Function
    T = t2 - t1
    ω₀ = 2π / T
    n = Int((2π * band) / ω₀)
    a₀ = 1 / T * quadgk(g, t1, t2)[1]
    an = zeros(Float64, n)
    bn = zeros(Float64, n)

    for i in 1:n
        an[i] = 2 / T * quadgk(t -> g(t) * cos(ω₀ * i * t), t1, t2)[1]
        bn[i] = 2 / T * quadgk(t -> g(t) * sin(ω₀ * i * t), t1, t2)[1]
    end

    function fourier_series(t)
        result = a₀ / 2
        for i in 1:n
            result += an[i] * cos(ω₀ * i * t) + bn[i] * sin(ω₀ * i * t)
        end
        return result
    end

    return fourier_series
end


function rand_signal_bl(f1::Real, f2::Real)::Function
    f = f1.+(f2-f1)*rand(1000)
    ϕ = -π .+ 2π .* rand(1000)
    A = randn(1000)
    A ./= sqrt(sum(A .^ 2) / 1000)
    
    return t -> sum(@. A * sin(2π * f * t + ϕ))    
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
    N = -M÷2:M÷2
    y = zeros(Float64,length(x))
    for n in 1:length(x)
        for m in N
            if n+m > 0 && n+m <=length(x)
            y[n]+=x[n+m]
            end
        end
    end
    return y/M
end

function running_energy(x::AbstractVector, M::Integer)::Vector
    N = -M÷2:M÷2
    y = zeros(Float64,length(x))
    for n in 1:length(x)
        for m in N
            if n+m > 0 && n+m <=length(x)
            y[n]+=x[n+m]^2
            end
        end
    end
    return y
end

function running_power(x::AbstractVector, M::Integer)::Vector
    N = -M÷2:M÷2
    y = zeros(Float64,length(x))
    for n in 1:length(x)
        for m in N
            if n+m > 0 && n+m <=length(x)
            y[n]+=x[n+m]^2
            end
        end
    end
    return y/M
end

# Próbkowanie
function interpolate(
    m::AbstractVector,
    s::AbstractVector,
    kernel::Function=sinc
)
    return t -> begin
        sum = 0
        Δt = m[2]-m[1]
        for n in eachindex(s)
            sum += kernel((t-m[n])/Δt)*s[n]
        end
        return sum
    end
end

# Kwantyzacja
quantize(L::AbstractVector)::Function = x -> L[argmin(abs.(-L .+ x))]
SQNR(N::Integer)::Real = 1.76 + 6.02*N
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
   N = length(x)
   w = OffsetArray(
        [exp(-1im * 2 * π/N * n) for n in 0:N-1],
        0:N-1
   )

   [
    sum((
        x[n+1] * w[(n*k)%N] for n in 0:N-1
        )) for k in 0:(N÷2)
   ]

end

function irdft(X::AbstractVector, N::Integer)::Vector
    S = length(X)
    X = [n <= S ? X[n] : conj(X[2S-n+(N % 2 == 0 ? 0 : 1)]) for n in 1:N]
    real.(idft(X))
end

function fft_radix2_dit_r(x::AbstractVector)::Vector
    function bitreverse(n, bits)
        reversed = 0
        for i in 1:bits
            reversed <<= 1
            reversed |= (n & 1)
            n >>= 1
        end
        return reversed
    end

    N = length(x)
    bits = Int(log2(N))
    X = copy(ComplexF64.(x))

    # Bit-reversal permutation
    for i in 1:N
        j = bitreverse(i - 1, bits) + 1
        if i < j
            X[i], X[j] = X[j], X[i]
        end
    end

    n₁ = 0
    n₂ = 1
    for i=1:bits
      n₁ = n₂ 
      n₂ *= 2
      
      step_angle = -2π/n₂
      angle = 0
      for j=1:n₁
        factors = exp(im*angle)
        angle += step_angle
        
        for k=j:n₂:N
          X[k], X[k+n₁] = X[k] + factors * X[k+n₁], X[k] - factors * X[k+n₁]
        end
      end
    end
    
    return X   
end

function ifft_radix2_dif_r(X::AbstractVector)::Vector
    X_R = [Complex(imag(z), real(z)) for z in X]
    x_r = fft_radix2_dit_r(X_R)
    x_n = [Complex(imag(z), real(z)) for z in x_r]
    return x_n/length(x_n)
end


# this is absolutely not optimised, but works
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
    X_R = [Complex(imag(z), real(z)) for z in X]
    x_r = fft(X_R)
    x_n = [Complex(imag(z), real(z)) for z in x_r]
    return x_n/length(x_n)
end


fftfreq(N::Int, fp::Float64) = fp * vcat(0:(N ÷ 2), -((N-1) ÷ 2):-1) / N
rfftfreq(N::Integer, fs::Real)::Vector = (0:N÷2) * fs / N
amplitude_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = abs.(fft(x.*w))/(length(x)*mean(w))
power_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = abs.(fft(x.*w).^2)/(length(x)*energy(w))
psd(x::AbstractVector, w::AbstractVector=rect(length(x)), fs::Real=1.0)::Vector = abs.(fft(x.*2).^2)/(length(x)*energy(w)*fs)

function periodogram(
    x::AbstractVector,
    w::AbstractVector=rect(length(x)),
    L::Integer = 0,
    fs::Real=1.0)::Vector
    N = length(x)
    K = length(w)
    step = K - L
    num_segments = div(N - L, step) + 1
    Pxx = zeros(Float64, K)

    for i in 0:(num_segments-1)
        start = i * step + 1
        segment = x[start:min(start+K-1,N)]

        if length(segment) < K
            segment = vcat(segment, zeros(K - length(segment)))
        end

        segment_psd = psd(segment, w, fs)

        Pxx += segment_psd

    end
    
    return Pxx / num_segments
end



function stft(x::AbstractVector, w::AbstractVector, L::Integer)::Matrix
    N = length(x)
    K = length(w)
    step = K - L
    num_segments = div(N - L, step) + 1
    stft_matrix = zeros(ComplexF64,(K ÷ 2 + 1,num_segments))
    for i in 0:(num_segments-1)
        start = i * step + 1
        segment = x[start:min(start+K-1,N)]
        if length(segment) < K
            segment = vcat(segment, zeros(K - length(segment)))
        end

        segment = segment .* w
        segment_fft = rfft(segment)
        stft_matrix[:,i+1] = segment_fft
    end
    return stft_matrix
end


function istft(X::AbstractMatrix{<:Complex}, w::AbstractVector{<:Real}, L::Integer)::AbstractVector{<:Real}
    K = length(w)
    step = K - L
    num_segments = size(X, 2)
    N = (num_segments - 1) * step + K

    x_reconstructed = zeros(Float64, N)
    window_sum = zeros(Float64, N)

    for i in 0:(num_segments-1)
        start = i * step + 1
        segment_ifft = irfft(X[:,i+1], K)

        x_reconstructed[start:start+K-1] += segment_ifft .* w
        window_sum[start:start+K-1] += w
    end

    x_reconstructed ./= window_sum
    
    return x_reconstructed
end

function conv(f::Vector, g::Vector)::Vector
    N = length(f)
    M = length(g)
    K = N + M - 1
    y = zeros(K)
    for n = 1:K
        for m = 1:M
            if(n-m+1 > 0) && (n-m+1 <= M)
                y[n] += f[m]*g[n-m+1]
            end
        end
    end

    return y
end

function fast_conv(f::Vector, g::Vector)::Vector
    N = length(f) + length(g) - 1

    f_padded = vcat(f, zeros(N - length(f)))
    g_padded = vcat(g, zeros(N - length(g)))

    F = fft(f_padded)
    G = fft(g_padded)
    Y = F .* G
    y = real(ifft(Y))

    return y
end

function overlap_add(x::Vector, h::Vector, L::Integer)::Vector
    N = length(x)
    M = length(h)
    P = L + M - 1
    
    result = zeros(Float64, N + M - 1)
    
    for i in 1:L:N
        segment = x[i:min(i+L-1, N)]
        
        conv_segment = fast_conv(segment, h)
        
        start_idx = i
        end_idx = min(i + P - 1, length(result))
        result[start_idx:end_idx] += conv_segment[1:(end_idx - start_idx + 1)]
    end
    
    return result
end

function overlap_save(x::Vector, h::Vector, L::Integer)::Vector
    M = length(h)
    N = L + M - 1

    padded_h = vcat(h, zeros(N - M))
    H = fft(padded_h)

    y = []
    padded_x = vcat(zeros(M - 1), x, zeros(N - 1))

    for k in 1:L:(length(padded_x)-N+1)
        xk = padded_x[k:k+N-1]
        Xk = fft(xk)
        Yk = ifft(H .* Xk)
        y = vcat(y, real(Yk[M:end]))
    end

    return y[1:(length(x)+M-1)]
end

function lti_filter(b::Vector, a::Vector, x::Vector)::Vector
    N = length(x)
    K = length(a)
    a /= a[1]
    M = length(b)
    y = zeros(N)
    for n in 1:N
        for m in 1:M
            if n-m+1 > 0 && n-m+1<=N
                y[n]+=b[m]*x[n-m+1]
            end
        end
        for k in 2:K
            if n-k+1>0 && n-k+1<=N
                y[n]-=a[k]*y[n-k+1]
            end
        end
    end
    return y
end

function filtfilt(b::Vector, a::Vector, x::Vector)::Vector
    y = zeros(Float64,length(x))
    for n in 1:length(x)
        for m in 1:length(b)
            if n-m+1 > 0 && n-m+1 <= length(x)
                y[n]+=b[m]*x[n-m+1]
            end
        end
        for k in 2:length(a)
            if n-k+1>0 && n-k+1<= length(y)
                y[n]-=a[k]*y[n-k+1]
            end
        end
    end

    return y
end

function lti_amp(f::Real, b::Vector, a::Vector)::Real
    licznik = 0
    mianownik = 0
    for l in 1:length(b)
        licznik+=b[l]*cispi(2f)^-(l-1)
    end
    for m in 1:length(a)
        mianownik+=a[m]*cispi(2f)^-(m-1)
    end
    return abs(licznik/mianownik)
end

function lti_phase(f::Real, b::Vector, a::Vector)::Real
    licznik = 0
    mianownik = 0
    for l in 1:length(b)
        licznik+=b[l]*cispi(2f)^-(l-1)
    end
    for m in 1:length(a)
        mianownik+=a[m]*cispi(2f)^-(m-1)
    end
    return angle(licznik/mianownik)
end


function firwin_lp_I(order, F0)
    n = -order/2:order/2
    h = [2F0*sinc(2*F0*n[i]) for i in eachindex(n)]
    return h
end

function firwin_hp_I(order, F0)
    n = -order/2:order/2
    h = [kronecker(n[i]) - 2F0*sinc(2F0*n[i]) for i in eachindex(n)]
    return h
end

function firwin_bp_I(order, F1, F2)
    n = -order/2:order/2
    h = [2F2*sinc(2F2*n[i]) - 2F1*sinc(2F1*n[i]) for i in eachindex(n)]
    return h
end

function firwin_bs_I(order, F1, F2)
    n = -order/2:order/2
    h = [kronecker(n[i]) - (2F2*sinc(2F2*n[i]) - 2F1*sinc(2F1*n[i])) for i in eachindex(n)]
    return h
end

function firwin_lp_II(N, F0)
    missing
end

function firwin_bp_II(N, F1, F2)
    missing
end

function firwin_diff(N::Int)
    x = -N:1:N
    h = [cospi(n)/n for n in x]
    return h
end

function resample(x::Vector, M::Int, N::Int)
    missing
end

end
