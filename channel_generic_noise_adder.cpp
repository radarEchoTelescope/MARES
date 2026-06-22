// =============================================================================
//  channel_generic_noise_adder.cpp
//
//  Standalone C++17 re-implementation of NuRadioReco's
//  `channelGenericNoiseAdder` Python module, using FFTW3 for the real<->
//  complex transforms (mirrors NumPy's `rfft`/`irfft`, which the original
//  Python relies on under the hood).
//
//  The module generates band-limited noise traces (not based on measured
//  data) that can be added to detector channel waveforms. Two spectral
//  shapes are supported:
//      - "perfect_white": flat amplitude spectrum
//      - "rayleigh":       per-bin amplitude drawn from a Rayleigh distribution
//  Random phases are assigned to every frequency bin (except DC, and the
//  Nyquist bin for even-length traces, which must stay real), and the
//  resulting half-spectrum is inverse-(real-)FFT'd back to the time domain.
//
//  Build:
//      g++ -std=c++17 -O2 channel_generic_noise_adder.cpp -lfftw3 -o noise_adder
//  (Debian/Ubuntu: `sudo apt-get install libfftw3-dev` to get the headers/lib.)
//
//  Design notes / faithfulness to the original Python:
//  -----------------------------------------------------------------------
//  * `fftutil::rfft`/`irfft` wrap `fftw_plan_dft_r2c_1d` / `fftw_plan_dft_c2r_1d`,
//    matching NumPy's `np.fft.rfft` / `np.fft.irfft` convention (forward
//    transform unnormalized, inverse transform divided by N). This mirrors
//    `NuRadioReco.utilities.fft`'s `freqs` / `time2freq` / `freq2time` helpers.
//  * The original Python uses NumPy's `Generator(Philox(seed))` counter-based
//    bit generator for reproducible, parallel-safe random streams. Exactly
//    reproducing Philox4x32-10 bit-for-bit is out of scope for a "standalone"
//    translation, so this version uses `std::mt19937_64` instead. The
//    *algorithm* (uniform phase draws, inverse-CDF Rayleigh sampling) is
//    identical; only the underlying bitstream differs, so output will be
//    statistically but not bit-for-bit identical to the Python version.
//  * Two latent bugs in the original Python were found while translating
//    and are fixed here (each is flagged at the call site below):
//      1. `bandlimited_noise_from_precalculated_parameters` reads
//         `self.precalculated_parameters["amplitude"]` for the
//         'perfect_white' branch, but `precalculate_bandlimited_noise_parameters`
//         never stores that key, so the original would raise a KeyError.
//         Fixed by actually storing `amplitude` in the precalculated struct.
//      2. `bandlimited_noise_from_spectrum`'s 'perfect_white' branch does
//         `ampl = amplitude * sigscale`, which rebinds the local Python name
//         to a bare scalar instead of filling the `ampl` array (so phases
//         could never be applied). Fixed by filling the whole array, which
//         matches the evident intent (consistent with the `bandlimited_noise`
//         method's perfect_white branch).
// =============================================================================

#include <fftw3.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <functional>
#include <iostream>
#include <map>
#include <optional>
#include <random>
#include <stdexcept>
#include <vector>

// -----------------------------------------------------------------------
//  units
//
//  A tiny self-consistent unit system (mirrors NuRadioReco.utilities.units
//  in spirit: every physical quantity is stored "in its base unit", and
//  named constants let you write expressive code like `50 * units::MHz`).
// -----------------------------------------------------------------------
namespace units {
constexpr double V = 1.0;
constexpr double mV = 1e-3 * V;

constexpr double Hz = 1.0;
constexpr double kHz = 1e3 * Hz;
constexpr double MHz = 1e6 * Hz;
constexpr double GHz = 1e9 * Hz;
}  // namespace units

// -----------------------------------------------------------------------
//  fftutil: FFTW3-backed, NuRadioReco-style real-FFT helpers
//  (freqs / time2freq / freq2time), matching NumPy's rfft/irfft convention.
// -----------------------------------------------------------------------
namespace fftutil {

using cplx = std::complex<double>;

// numpy.fft.rfftfreq(n, 1/sampling_rate) equivalent: returns the n/2+1
// non-negative frequency bins of a real signal of length n sampled at
// `sampling_rate`.
inline std::vector<double> freqs(int n_samples, double sampling_rate) {
    const int half_len = n_samples / 2 + 1;
    std::vector<double> f(half_len);
    for (int k = 0; k < half_len; ++k) {
        f[k] = static_cast<double>(k) * sampling_rate / static_cast<double>(n_samples);
    }
    return f;
}

// numpy.fft.rfft equivalent: forward real-to-complex FFT, n/2+1 complex
// bins, unnormalized (X[k] = sum_n x[n] exp(-2*pi*i*k*n/N)).
inline std::vector<cplx> rfft(const std::vector<double>& trace) {
    const int n = static_cast<int>(trace.size());
    const int half_len = n / 2 + 1;

    std::vector<double> in(trace);  // local scratch copy; r2c may overwrite the input buffer
    std::vector<fftw_complex> out(half_len);

    fftw_plan plan = fftw_plan_dft_r2c_1d(n, in.data(), out.data(), FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    std::vector<cplx> result(half_len);
    for (int k = 0; k < half_len; ++k) result[k] = cplx(out[k][0], out[k][1]);
    return result;
}

// numpy.fft.irfft equivalent: complex-to-real inverse FFT from the n/2+1
// supplied half-spectrum bins (Hermitian symmetry is implicit, exactly as
// FFTW's c2r transform assumes), normalized by 1/n, returns length `n`.
inline std::vector<double> irfft(const std::vector<cplx>& half_spectrum, int n) {
    const int half_len = n / 2 + 1;

    std::vector<fftw_complex> in(half_len);
    for (int k = 0; k < half_len; ++k) {
        if (k < static_cast<int>(half_spectrum.size())) {
            in[k][0] = half_spectrum[k].real();
            in[k][1] = half_spectrum[k].imag();
        } else {
            in[k][0] = 0.0;
            in[k][1] = 0.0;
        }
    }
    std::vector<double> out(n);

    // FFTW's c2r transform destroys its input by default; that's fine since
    // `in` here is a local scratch copy.
    fftw_plan plan = fftw_plan_dft_c2r_1d(n, in.data(), out.data(), FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (double& v : out) v /= static_cast<double>(n);  // FFTW inverse transforms are unnormalized
    return out;
}

// NuRadioReco convention: time2freq divides by sampling_rate, freq2time
// multiplies by it, so that spectral amplitudes carry consistent units
// (e.g. V / GHz) across the transform.
inline std::vector<cplx> time2freq(const std::vector<double>& trace, double sampling_rate) {
    auto spec = rfft(trace);
    for (auto& v : spec) v /= sampling_rate;
    return spec;
}

inline std::vector<double> freq2time(const std::vector<cplx>& spectrum, double sampling_rate, int n) {
    auto trace = irfft(spectrum, n);
    for (auto& v : trace) v *= sampling_rate;
    return trace;
}

// scipy.integrate.trapezoid(y, x) equivalent.
inline double trapz(const std::vector<double>& y, const std::vector<double>& x) {
    double sum = 0.0;
    for (size_t i = 1; i < y.size(); ++i) {
        sum += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
    }
    return sum;
}

}  // namespace fftutil

// -----------------------------------------------------------------------
//  channelGenericNoiseAdder
// -----------------------------------------------------------------------
enum class NoiseType { PerfectWhite, Rayleigh };

class ChannelGenericNoiseAdder {
public:
    // Holds the state cached by PrecalculateBandlimitedNoiseParameters(),
    // mirroring the Python `self.precalculated_parameters` dict.
    struct PrecalculatedParams {
        int n_samples_freq = 0;
        std::vector<bool> selection;
        int nbinsactive = 0;
        double sigscale = 0.0;
        double fsigma = 0.0;
        double amplitude = 0.0;  // NOTE: bug fix vs. original Python, see file header (#1)
        double sampling_rate = 0.0;
        std::vector<double> frequencies;
        int n_samples = 0;
    };

    ChannelGenericNoiseAdder() { Begin(); }

    // ---------------------------------------------------------------
    // begin()
    // ---------------------------------------------------------------
    void Begin(bool debug = false, std::optional<uint64_t> seed = std::nullopt) {
        debug_ = debug;
        if (seed.has_value()) {
            rng_.seed(*seed);
        } else {
            std::random_device rd;
            rng_.seed(rd());
        }
    }

    // ---------------------------------------------------------------
    // add_random_phases()
    //
    // Multiplies amps[1 .. Np] (inclusive) by a unit-modulus complex
    // factor with a uniformly random phase. Np = (n_samples_time_domain-1)/2
    // (integer division) is the number of strictly-positive frequency bins
    // of the *full* bilateral spectrum of length `n_samples_time_domain`;
    // this leaves the DC bin (and, for even n_samples_time_domain, the
    // Nyquist bin) untouched so the resulting half-spectrum can be fed
    // straight into an irfft.
    // ---------------------------------------------------------------
    std::vector<fftutil::cplx> AddRandomPhases(const std::vector<double>& amps, int n_samples_time_domain) {
        std::vector<fftutil::cplx> out(amps.size());
        for (size_t i = 0; i < amps.size(); ++i) out[i] = fftutil::cplx(amps[i], 0.0);
        return AddRandomPhasesInPlace(std::move(out), n_samples_time_domain);
    }

    std::vector<fftutil::cplx> AddRandomPhasesInPlace(std::vector<fftutil::cplx> amps, int n_samples_time_domain) {
        const int Np = (n_samples_time_domain - 1) / 2;
        std::uniform_real_distribution<double> uni(0.0, 2.0 * M_PI);
        const int upper = std::min<int>(Np, static_cast<int>(amps.size()) - 1);
        for (int i = 1; i <= upper; ++i) {
            double phase = uni(rng_);
            amps[i] *= fftutil::cplx(std::cos(phase), std::sin(phase));
        }
        return amps;
    }

    // ---------------------------------------------------------------
    // bandlimited_noise() — frequency-domain half-spectrum.
    //
    // `min_freq`/`max_freq`/`bandwidth` use std::optional to mirror
    // Python's `None` defaults.
    // ---------------------------------------------------------------
    std::vector<fftutil::cplx> BandlimitedNoiseFreqDomain(
        std::optional<double> min_freq_in, std::optional<double> max_freq_in, int n_samples, double sampling_rate,
        double amplitude, NoiseType type = NoiseType::PerfectWhite, std::optional<double> bandwidth = std::nullopt) {
        std::vector<double> frequencies = fftutil::freqs(n_samples, sampling_rate);
        const int n_samples_freq = static_cast<int>(frequencies.size());

        double min_freq;
        if (!min_freq_in.has_value() || *min_freq_in == 0.0) {
            // remove DC component; see Python docstring for the rationale
            // behind taking a bin-spacing difference rather than frequencies[1].
            min_freq = 0.5 * (frequencies[2] - frequencies[1]);
        } else {
            min_freq = *min_freq_in;
        }

        double max_freq;
        if (!max_freq_in.has_value()) {
            max_freq = frequencies.back();
        } else {
            max_freq = *max_freq_in;
            if (std::round(max_freq * 1000.0) / 1000.0 > std::round(frequencies.back() * 1000.0) / 1000.0) {
                std::cerr << "WARNING: max_freq (" << max_freq / units::MHz
                          << " MHz) is above the Nyquist frequency (" << frequencies.back() / units::MHz
                          << " MHz). The simulated noise amplitude might deviate from what you intended.\n";
            }
        }

        std::vector<bool> selection(n_samples_freq);
        int nbinsactive = 0;
        for (int i = 0; i < n_samples_freq; ++i) {
            selection[i] = (frequencies[i] >= min_freq) && (frequencies[i] <= max_freq);
            if (selection[i]) ++nbinsactive;
        }
        if (nbinsactive == 0) {
            throw std::runtime_error("BandlimitedNoiseFreqDomain: no frequency bins selected — check min/max_freq.");
        }

        if (bandwidth.has_value()) {
            double sampling_bandwidth = std::min(0.5 * sampling_rate, max_freq) - min_freq;
            amplitude *= 1.0 / std::sqrt(*bandwidth / sampling_bandwidth);
        }

        std::vector<double> ampl(n_samples_freq, 0.0);
        const double sigscale = static_cast<double>(n_samples) / std::sqrt(static_cast<double>(nbinsactive));

        if (type == NoiseType::PerfectWhite) {
            for (int i = 0; i < n_samples_freq; ++i) {
                if (selection[i]) ampl[i] = amplitude * sigscale;
            }
        } else {  // Rayleigh
            const double fsigma = amplitude * sigscale / std::sqrt(2.0);
            std::uniform_real_distribution<double> uni(0.0, 1.0);
            for (int i = 0; i < n_samples_freq; ++i) {
                if (selection[i]) ampl[i] = SampleRayleigh(fsigma, uni(rng_));
            }
        }

        auto noise = AddRandomPhases(ampl, n_samples);
        for (auto& v : noise) v /= sampling_rate;
        return noise;
    }

    // Time-domain convenience wrapper around BandlimitedNoiseFreqDomain().
    std::vector<double> BandlimitedNoiseTimeDomain(std::optional<double> min_freq, std::optional<double> max_freq,
                                                    int n_samples, double sampling_rate, double amplitude,
                                                    NoiseType type = NoiseType::PerfectWhite,
                                                    std::optional<double> bandwidth = std::nullopt) {
        auto noise_freq =
            BandlimitedNoiseFreqDomain(min_freq, max_freq, n_samples, sampling_rate, amplitude, type, bandwidth);
        return fftutil::freq2time(noise_freq, sampling_rate, n_samples);
    }

    // ---------------------------------------------------------------
    // precalculate_bandlimited_noise_parameters() /
    // bandlimited_noise_from_precalculated_parameters()
    //
    // Splits the (expensive, amplitude-independent) bin-selection work
    // from the (cheap, per-call) random draw, exactly like the Python
    // version — useful when generating many independent noise traces with
    // the same min/max frequency and sample count.
    // ---------------------------------------------------------------
    void PrecalculateBandlimitedNoiseParameters(std::optional<double> min_freq_in, std::optional<double> max_freq_in,
                                                 int n_samples, double sampling_rate, double amplitude,
                                                 std::optional<double> bandwidth = std::nullopt) {
        std::vector<double> frequencies = fftutil::freqs(n_samples, sampling_rate);
        const int n_samples_freq = static_cast<int>(frequencies.size());

        double min_freq;
        if (!min_freq_in.has_value() || *min_freq_in == 0.0) {
            min_freq = 0.5 * (frequencies[2] - frequencies[1]);
        } else {
            min_freq = *min_freq_in;
        }
        double max_freq = max_freq_in.has_value() ? *max_freq_in : frequencies.back();

        std::vector<bool> selection(n_samples_freq);
        int nbinsactive = 0;
        for (int i = 0; i < n_samples_freq; ++i) {
            selection[i] = (frequencies[i] >= min_freq) && (frequencies[i] <= max_freq);
            if (selection[i]) ++nbinsactive;
        }

        if (bandwidth.has_value()) {
            double sampling_bandwidth = std::min(0.5 * sampling_rate, max_freq) - min_freq;
            amplitude *= 1.0 / std::sqrt(*bandwidth / sampling_bandwidth);
        }

        const double sigscale = static_cast<double>(n_samples) / std::sqrt(static_cast<double>(nbinsactive));
        const double fsigma = amplitude * sigscale / std::sqrt(2.0);

        precalculated_ = PrecalculatedParams{n_samples_freq, selection,    nbinsactive, sigscale,
                                              fsigma,         amplitude,   sampling_rate, frequencies,
                                              n_samples};
    }

    std::vector<fftutil::cplx> BandlimitedNoiseFromPrecalculatedParametersFreqDomain(
        NoiseType type = NoiseType::PerfectWhite) {
        if (precalculated_.n_samples == 0) {
            throw std::runtime_error("Call PrecalculateBandlimitedNoiseParameters() first.");
        }
        const auto& p = precalculated_;
        std::vector<double> ampl(p.n_samples_freq, 0.0);

        if (type == NoiseType::PerfectWhite) {
            // amplitude is read from the cached struct here — see file
            // header note (#1): the original Python looks this up too,
            // but never stored it, so this branch would raise a KeyError
            // upstream. We store it in PrecalculateBandlimitedNoiseParameters().
            for (int i = 0; i < p.n_samples_freq; ++i) {
                if (p.selection[i]) ampl[i] = p.amplitude * p.sigscale;
            }
        } else {
            std::uniform_real_distribution<double> uni(0.0, 1.0);
            for (int i = 0; i < p.n_samples_freq; ++i) {
                if (p.selection[i]) ampl[i] = SampleRayleigh(p.fsigma, uni(rng_));
            }
        }

        auto noise = AddRandomPhases(ampl, p.n_samples);
        for (auto& v : noise) v /= p.sampling_rate;
        return noise;
    }

    std::vector<double> BandlimitedNoiseFromPrecalculatedParametersTimeDomain(
        NoiseType type = NoiseType::PerfectWhite) {
        auto noise_freq = BandlimitedNoiseFromPrecalculatedParametersFreqDomain(type);
        return fftutil::freq2time(noise_freq, precalculated_.sampling_rate, precalculated_.n_samples);
    }

    // ---------------------------------------------------------------
    // bandlimited_noise_from_spectrum()
    //
    // `spectrum` may be supplied directly as a vector of n/2+1 complex
    // bins, or as a function `double freq -> complex amplitude` evaluated
    // at each rfft frequency bin (mirroring Python's "array or callable"
    // duality).
    // ---------------------------------------------------------------
    std::vector<fftutil::cplx> BandlimitedNoiseFromSpectrumFreqDomain(
        int n_samples, double sampling_rate, std::vector<fftutil::cplx> spectrum,
        std::optional<double> amplitude_in = std::nullopt, NoiseType type = NoiseType::PerfectWhite) {
        std::vector<double> frequencies = fftutil::freqs(n_samples, sampling_rate);
        const int n_samples_freq = static_cast<int>(frequencies.size());

        std::vector<bool> selection(n_samples_freq);
        int n_active = 0;
        for (int i = 0; i < n_samples_freq; ++i) {
            selection[i] = frequencies[i] > 0.0;
            if (selection[i]) ++n_active;
        }

        double amplitude;
        double sigscale;
        if (amplitude_in.has_value()) {
            std::vector<double> spectrum_mag2(n_samples_freq);
            for (int i = 0; i < n_samples_freq; ++i) spectrum_mag2[i] = std::norm(spectrum[i]);
            double norm = fftutil::trapz(spectrum_mag2, frequencies);
            double max_freq = frequencies.back();
            amplitude = *amplitude_in / std::sqrt(norm / max_freq);
            sigscale = static_cast<double>(n_samples) / std::sqrt(static_cast<double>(n_active));
        } else {
            amplitude = std::sqrt(static_cast<double>(n_samples));
            sigscale = 1.0;
        }

        std::vector<double> ampl(n_samples_freq, 0.0);
        if (type == NoiseType::PerfectWhite) {
            // NOTE: bug fix vs. original Python, see file header (#2):
            // the source rebinds the local `ampl` name to a bare scalar
            // (`ampl = amplitude * sigscale`) instead of filling the array
            // (`ampl[:] = ...`), which would silently produce un-phased,
            // non-array noise. We fill the whole array, matching the
            // evident intent and the behaviour of bandlimited_noise()'s
            // own 'perfect_white' branch.
            std::fill(ampl.begin(), ampl.end(), amplitude * sigscale);
        } else {
            const double fsigma = amplitude * sigscale / std::sqrt(2.0);
            std::uniform_real_distribution<double> uni(0.0, 1.0);
            for (int i = 0; i < n_samples_freq; ++i) {
                if (selection[i]) ampl[i] = SampleRayleigh(fsigma, uni(rng_));
            }
        }

        auto noise = AddRandomPhases(ampl, n_samples);
        for (auto& v : noise) v /= sampling_rate;
        for (int i = 0; i < n_samples_freq; ++i) noise[i] *= spectrum[i];
        return noise;
    }

    std::vector<fftutil::cplx> BandlimitedNoiseFromSpectrumFreqDomain(
        int n_samples, double sampling_rate, const std::function<fftutil::cplx(double)>& spectrum_fn,
        std::optional<double> amplitude_in = std::nullopt, NoiseType type = NoiseType::PerfectWhite) {
        std::vector<double> frequencies = fftutil::freqs(n_samples, sampling_rate);
        std::vector<fftutil::cplx> spectrum(frequencies.size());
        for (size_t i = 0; i < frequencies.size(); ++i) spectrum[i] = spectrum_fn(frequencies[i]);
        return BandlimitedNoiseFromSpectrumFreqDomain(n_samples, sampling_rate, std::move(spectrum), amplitude_in,
                                                        type);
    }

    std::vector<double> BandlimitedNoiseFromSpectrumTimeDomain(
        int n_samples, double sampling_rate, std::vector<fftutil::cplx> spectrum,
        std::optional<double> amplitude_in = std::nullopt, NoiseType type = NoiseType::PerfectWhite) {
        auto noise_freq =
            BandlimitedNoiseFromSpectrumFreqDomain(n_samples, sampling_rate, std::move(spectrum), amplitude_in, type);
        return fftutil::freq2time(noise_freq, sampling_rate, n_samples);
    }

    // ---------------------------------------------------------------
    // run()
    //
    // Stand-ins for NuRadioReco's Event/Station/Detector: a Channel is
    // just an id + trace + sampling rate, a Station is a list of Channels.
    // ---------------------------------------------------------------
    struct Channel {
        int id;
        std::vector<double> trace;
        double sampling_rate;
    };
    struct Station {
        std::vector<Channel> channels;
    };

    void Run(Station& station, double amplitude = 1.0 * units::mV,
             std::optional<double> min_freq = 50.0 * units::MHz,
             std::optional<double> max_freq = 2000.0 * units::MHz, NoiseType type = NoiseType::PerfectWhite,
             const std::vector<int>& excluded_channels = {}, std::optional<double> bandwidth = std::nullopt,
             const std::optional<std::map<int, double>>& amplitude_per_channel = std::nullopt) {
        for (auto& channel : station.channels) {
            if (std::find(excluded_channels.begin(), excluded_channels.end(), channel.id) !=
                excluded_channels.end()) {
                continue;
            }

            double tmp_amplitude = amplitude;
            if (amplitude_per_channel.has_value()) {
                auto it = amplitude_per_channel->find(channel.id);
                if (it != amplitude_per_channel->end()) tmp_amplitude = it->second;
            }

            auto noise = BandlimitedNoiseTimeDomain(min_freq, max_freq, static_cast<int>(channel.trace.size()),
                                                     channel.sampling_rate, tmp_amplitude, type, bandwidth);

            if (debug_) {
                double rms = 0.0;
                for (double v : noise) rms += v * v;
                rms = std::sqrt(rms / noise.size());
                std::cerr << "DEBUG: input amplitude " << tmp_amplitude << ", noise RMS " << rms << "\n";
            }

            for (size_t i = 0; i < channel.trace.size(); ++i) channel.trace[i] += noise[i];
        }
    }

    void End() {}

private:
    // Inverse-CDF sampling for the Rayleigh distribution with scale `sigma`,
    // given a uniform(0,1) draw `u`.
    static double SampleRayleigh(double sigma, double u) {
        u = std::clamp(u, std::numeric_limits<double>::min(), 1.0 - std::numeric_limits<double>::epsilon());
        return sigma * std::sqrt(-2.0 * std::log(u));
    }

    bool debug_ = false;
    std::mt19937_64 rng_;
    PrecalculatedParams precalculated_;
};

// =============================================================================
//  Demonstration / smoke test
// =============================================================================
#ifndef CHANNEL_GENERIC_NOISE_ADDER_NO_MAIN
int main() {
    ChannelGenericNoiseAdder noise_adder;
    noise_adder.Begin(/*debug=*/false, /*seed=*/12345);

    const int n_samples = 1024;
    const double sampling_rate = 3.2 * units::GHz;

    // --- perfect_white noise, generated directly ---------------------------
    auto white_noise = noise_adder.BandlimitedNoiseTimeDomain(
        /*min_freq=*/50 * units::MHz, /*max_freq=*/1000 * units::MHz, n_samples, sampling_rate,
        /*amplitude=*/1 * units::mV, NoiseType::PerfectWhite);

    double rms = 0.0;
    for (double v : white_noise) rms += v * v;
    rms = std::sqrt(rms / white_noise.size());
    std::cout << "perfect_white noise: " << n_samples << " samples, RMS = " << rms / units::mV << " mV "
              << "(target ~1 mV)\n";

    // --- rayleigh noise, generated directly ---------------------------------
    auto rayleigh_noise = noise_adder.BandlimitedNoiseTimeDomain(50 * units::MHz, 1000 * units::MHz, n_samples,
                                                                  sampling_rate, 1 * units::mV, NoiseType::Rayleigh);
    rms = 0.0;
    for (double v : rayleigh_noise) rms += v * v;
    rms = std::sqrt(rms / rayleigh_noise.size());
    std::cout << "rayleigh noise:      " << n_samples << " samples, RMS = " << rms / units::mV << " mV "
              << "(target ~1 mV)\n";

    // --- precalculated-parameters workflow ----------------------------------
    noise_adder.PrecalculateBandlimitedNoiseParameters(50 * units::MHz, 1000 * units::MHz, n_samples, sampling_rate,
                                                        1 * units::mV);
    auto precalc_noise = noise_adder.BandlimitedNoiseFromPrecalculatedParametersTimeDomain(NoiseType::PerfectWhite);
    rms = 0.0;
    for (double v : precalc_noise) rms += v * v;
    rms = std::sqrt(rms / precalc_noise.size());
    std::cout << "precalculated noise: " << n_samples << " samples, RMS = " << rms / units::mV << " mV "
              << "(target ~1 mV)\n";

    // --- run() over a small station ------------------------------------------
    ChannelGenericNoiseAdder::Station station;
    station.channels.push_back({0, std::vector<double>(n_samples, 0.0), sampling_rate});
    station.channels.push_back({1, std::vector<double>(n_samples, 0.0), sampling_rate});

    std::map<int, double> per_channel_amplitude{{0, 1 * units::mV}, {1, 2 * units::mV}};
    noise_adder.Run(station, 1 * units::mV, 50 * units::MHz, 1000 * units::MHz, NoiseType::PerfectWhite,
                     /*excluded_channels=*/{}, /*bandwidth=*/std::nullopt, per_channel_amplitude);

    for (auto& channel : station.channels) {
        double channel_rms = 0.0;
        for (double v : channel.trace) channel_rms += v * v;
        channel_rms = std::sqrt(channel_rms / channel.trace.size());
        std::cout << "station channel " << channel.id << ": RMS = " << channel_rms / units::mV << " mV\n";
    }

    noise_adder.End();
    return 0;
}
#endif
