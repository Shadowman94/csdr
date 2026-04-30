/*
This software is part of libcsdr, a set of simple DSP routines for
Software Defined Radio.

Copyright (c) 2014, Andras Retzler <randras@sdr.hu>
Copyright (c) 2019-2021 Jakob Ketterl <jakob.ketterl@gmx.de>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the copyright holder nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ANDRAS RETZLER BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#include "module.hpp"
#include "complex.hpp"
#include "fir.hpp"
#include "fmdemod.hpp"

#include <queue>
#include <functional>
#include <future>
#include <iostream>

namespace Csdr {

    class PilotPLL {
    private:
        double samplerate_;  // Sampling rate in Hz
        double pilot_freq_;  // Pilot frequency (typically 19 kHz)
        double alpha_;       // Integral gain (wn^2)
        double beta_;        // Proportional gain (2*damp*wn)
        double minfreq_;     // Minimum frequency (phase inc)
        double maxfreq_;     // Maximum frequency (phase inc)
        double phzref_;      // Reference phase increment
        double freq_;        // Current frequency (phase inc)
        double phase_;       // Current phase of the NCO
        
        // Lock detection (from original)
        double lock_;        // Lock metric
        double lockalpha_;   // Lock filter alpha
        double lockbeta_;    // Lock filter beta
        double locklimit_;   // Lock threshold
        double lockdelay_;   // Lock delay time in samples
        double lockcount_;   // Lock delay counter
        
    public:
        // Constructor
        // sample_rate: Audio sampling rate (e.g., 44100 Hz)
        // pilot_freq: Pilot tone frequency (default 19000 Hz)
        // damp: Damping factor (default 0.707 for critical damping)
        // bw: Loop bandwidth (e.g., 20 Hz for stable lock)
        PilotPLL(double sample_rate, double pilot_freq = 19000.0, double damp = 0.707, double bw = 30.0)
        : samplerate_(sample_rate), pilot_freq_(pilot_freq), phase_(0.0), lock_(0.0), lockcount_(0.0) {
            
            double fn = bw * 0.707;  // Natural frequency approximation
            double wn = 2.0 * M_PI * fn / samplerate_;
            
            alpha_ = wn * wn;            // Integral gain
            beta_ = 2.0 * damp * wn;     // Proportional gain
            
            double phz = 2.0 * M_PI * pilot_freq / samplerate_;
            // SDR++-style pilot lock window: 18.75..19.25 kHz (±250 Hz).
            minfreq_ = phz - 2.0 * M_PI * 250.0 / samplerate_;
            maxfreq_ = phz + 2.0 * M_PI * 250.0 / samplerate_;
            phzref_ = phz;
            freq_ = phz;
            
            // Lock detection parameters
            lockalpha_ = 1.0 - std::exp(-1.0 / (samplerate_ * 0.2));  // 0.2 sec time constant
            lockbeta_ = 1.0 - lockalpha_;
            locklimit_ = 0.1;  // Lock threshold (low error means locked)
            lockdelay_ = samplerate_ * 0.5;  // 0.5 sec delay for stability
        }
        
        // Complex (I/Q) pilot PLL. This mirrors SDR++ behavior better than
        // the legacy real-only detector and is significantly less prone to
        // polarity ambiguity under jitter/load.
        double processIQ(double input_i, double input_q, double& pilot_strength) {
            double nco_i = std::cos(phase_);
            double nco_q = std::sin(phase_);

            // error = angle(input * conj(nco))
            double re = input_i * nco_i + input_q * nco_q;
            double im = input_q * nco_i - input_i * nco_q;
            double error = std::atan2(im, re);

            if (error > 0.7) error = 0.7;
            else if (error < -0.7) error = -0.7;

            freq_ += alpha_ * error;
            phase_ += freq_ + beta_ * error;

            if (phase_ > 2.0 * M_PI) phase_ -= 2.0 * M_PI;
            else if (phase_ < -2.0 * M_PI) phase_ += 2.0 * M_PI;

            if (freq_ > maxfreq_) freq_ = maxfreq_;
            else if (freq_ < minfreq_) freq_ = minfreq_;

            double mag = std::sqrt(input_i * input_i + input_q * input_q);
            lock_ = lock_ * lockalpha_ + lockbeta_ * std::fabs(error);
            if (lock_ < locklimit_ && mag > 1e-4) {
                lockcount_ = lockdelay_;
            }
            if (lockcount_ > 0.0) lockcount_--;

            pilot_strength = (lockcount_ > 0.0) ? (1.0 - lock_) : 0.0;
            if (pilot_strength < 0.0) pilot_strength = 0.0;
            if (pilot_strength > 1.0) pilot_strength = 1.0;

            return std::cos(2.0 * phase_);
        }

        double process(double input, double& pilot_strength) {
            return processIQ(input, 0.0, pilot_strength);
        }
        
        void reset() {
            phase_ = 0.0;
            freq_ = phzref_;
            lock_ = 0.0;
            lockcount_ = 0.0;
        }
        
        double getPhase() {
            return phase_;
        }
    };

    class BiquadFilter {
    private:
        double b0, b1, b2, a1, a2;  // Filter coefficients
        double x1, x2, y1, y2;     // Delay line states
        
    public:
        BiquadFilter() : b0(1.0), b1(0.0), b2(0.0), a1(0.0), a2(0.0), 
                        x1(0.0), x2(0.0), y1(0.0), y2(0.0) {}
        
        void setBandpass(double fc, double Q, double fs) {
            double omega = 2.0 * M_PI * fc / fs;
            double alpha = sin(omega) / (2.0 * Q);
            double cos_omega = cos(omega);
            
            // Bandpass coefficients
            double norm = 1.0 + alpha;
            b0 = alpha / norm;
            b1 = 0.0;
            b2 = -alpha / norm;
            a1 = -2.0 * cos_omega / norm;
            a2 = (1.0 - alpha) / norm;
        }

        void setBandpass2(double fc, double bw, double fs) {
            double Q = fc / bw;
            double omega = 2.0 * M_PI * fc / fs;
            // Optional: Adjust for bilinear warping if fc is high relative to fs/2
            double adj = (omega > 0.0) ? std::sin(omega) / omega : 1.0;
            Q /= adj;  // Compensates bandwidth cramping; omit if fs >> 2*fc (e.g., 192 kHz+)
            double alpha = std::sin(omega) / (2.0 * Q);
            double cos_omega = std::cos(omega);

            double norm = 1.0 + alpha;
            b0 = alpha / norm;
            b1 = 0.0;
            b2 = -alpha / norm;
            a1 = -2.0 * cos_omega / norm;
            a2 = (1.0 - alpha) / norm;
        }
        
        void setLowpassWithQ(double fc, double Q, double fs) {
            double omega = 2.0 * M_PI * fc / fs;
            double cos_omega = cos(omega);
            double alpha = sin(omega) / (2.0 * Q);
            
            double norm = 1.0 + alpha;
            b0 = (1.0 - cos_omega) / 2.0 / norm;
            b1 = (1.0 - cos_omega) / norm;
            b2 = (1.0 - cos_omega) / 2.0 / norm;
            a1 = -2.0 * cos_omega / norm;
            a2 = (1.0 - alpha) / norm;
        }
        
        void setLowpass(double fc, double fs) {
            setLowpassWithQ(fc, 0.707, fs);  // Standard Butterworth Q
        }
        
        void setHighpass(double fc, double fs) {
            double omega = 2.0 * M_PI * fc / fs;
            double cos_omega = cos(omega);
            double alpha = sin(omega) / sqrt(2.0);
            
            double norm = 1.0 + alpha;
            b0 = (1.0 + cos_omega) / 2.0 / norm;
            b1 = -(1.0 + cos_omega) / norm;
            b2 = (1.0 + cos_omega) / 2.0 / norm;
            a1 = -2.0 * cos_omega / norm;
            a2 = (1.0 - alpha) / norm;
        }
        
        double process(double input) {
            double output = b0 * input + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
            
            // Update delay line
            x2 = x1;
            x1 = input;
            y2 = y1;
            y1 = output;
            
            return output;
        }
        
        void reset() {
            x1 = x2 = y1 = y2 = 0.0;
        }
    };

    class MultistageFilter {
    private:
        // 12 stages = up to 24th-order Butterworth (more than enough headroom).
        // For 16th-order the highest required Q is ≈ 5.1 — well inside the clamp below.
        static const int MAX_STAGES = 12;
        BiquadFilter stages[MAX_STAGES];
        int num_stages;
        
    public:
        MultistageFilter() : num_stages(1) {}
        
        void setLowpass(double fc, double fs, int order = 6) {
            num_stages = (order + 1) / 2;  // Each biquad = 2 poles
            if (num_stages > MAX_STAGES) num_stages = MAX_STAGES;
            
            // Calculate Butterworth poles for higher order filter
            for (int i = 0; i < num_stages; ++i) {
                double angle = M_PI * (2 * i + 1) / (2 * order);
                double Q = 1.0 / (2.0 * cos(angle));
                
                // Limit Q to reasonable values
                if (Q < 0.5) Q = 0.5;
                if (Q > 10.0) Q = 10.0;
                
                stages[i].setLowpassWithQ(fc, Q, fs);
            }
        }
        
        void setBandpass(double fc, double Q, double fs) {
            num_stages = 1;
            stages[0].setBandpass(fc, Q, fs);
        }
        
        void setHighpass(double fc, double fs) {
            num_stages = 1;
            stages[0].setHighpass(fc, fs);
        }
        
        double process(double input) {
            double output = input;
            for (int i = 0; i < num_stages; ++i) {
                output = stages[i].process(output);
            }
            return output;
        }
        
        void reset() {
            for (int i = 0; i < num_stages; ++i) {
                stages[i].reset();
            }
        }
    };
    
    class NotchFilter {
    private:
        double b0, b1, b2, a1, a2;
        double x1, x2, y1, y2;
        
    public:
        NotchFilter() : b0(1.0), b1(0.0), b2(0.0), a1(0.0), a2(0.0),
                        x1(0.0), x2(0.0), y1(0.0), y2(0.0) {}
        
        void setNotch(double fc, double Q, double fs) {
            double omega = 2.0 * M_PI * fc / fs;
            double alpha = sin(omega) / (2.0 * Q);
            double cos_omega = cos(omega);
            
            double norm = 1.0 + alpha;
            b0 = 1.0 / norm;
            b1 = -2.0 * cos_omega / norm;
            b2 = 1.0 / norm;
            a1 = -2.0 * cos_omega / norm;
            a2 = (1.0 - alpha) / norm;
        }
        
        double process(double input) {
            double output = b0 * input + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
            
            x2 = x1; x1 = input;
            y2 = y1; y1 = output;
            
            return output;
        }
        
        void reset() {
            x1 = x2 = y1 = y2 = 0.0;
        }
    };
        
    class SimpleDCBlock {
        float R;
        float prev_in;
        float prev_out;
    public:
        SimpleDCBlock(float R_ = 0.995f) : R(R_), prev_in(0.0f), prev_out(0.0f) {}
        float process(float x) {
            float y = x - prev_in + R * prev_out;
            prev_in = x;
            prev_out = y;
            return y;
        }
    };

    template <typename T>
    class MonoFractionalDecimator {

        public:
            MonoFractionalDecimator();
            ~MonoFractionalDecimator();
            
            struct Denominator {

                Denominator() { 
                    readerAvailable = 0;
                    writterAvailable = 0;

                    where = 0;
                    num_poly_points = 0;

                    poly_precalc_denomiator.resize(0);
                    coeffs_buf.resize(0);

                    xifirst = 0;
                    xilast = 0;
                    rate = 0;
                    output_processed = 0;
                }

                size_t readerAvailable;
                size_t writterAvailable;

                double where;
                size_t output_processed;
                unsigned int num_poly_points; //number of samples that the Lagrange interpolator will use
                std::vector<float> poly_precalc_denomiator; //while we don't precalculate coefficients here as in a Farrow structure, because it is a fractional interpolator, but we rather precaculate part of the interpolator expression
                std::vector<float> coeffs_buf;
                int xifirst;
                int xilast;
                double rate;
            };

            struct DenominatorImmutable {
                unsigned int num_poly_points;
                std::vector<float> poly_precalc_denomiator;
                int xifirst;
                int xilast;
                double rate;

                FirFilter<T, float>* filter;

                DenominatorImmutable(unsigned int n, double r, FirFilter<T, float>* f) 
                    : num_poly_points(n & ~1), rate(r), filter(f)
                {
                    xifirst = -(num_poly_points / 2) + 1;
                    xilast  = num_poly_points / 2;
                    poly_precalc_denomiator.resize(num_poly_points);

                    for (int xi = xifirst, id = 0; xi <= xilast; xi++, id++) {
                        poly_precalc_denomiator[id] = 1;
                        for (int xj = xifirst; xj <= xilast; xj++) {
                            if (xi != xj) poly_precalc_denomiator[id] *= (xi - xj);
                        }
                    }
                }
            };

            struct DenominatorState {
                double where; // szál-specifikus állapot
                size_t output_processed;
                std::vector<float> coeffs_buf;

                std::mutex processMutex;

                DenominatorState() {};

                DenominatorState(int startWhere, int num_poly_points) 
                    : where(startWhere), output_processed(0) {
                        coeffs_buf.resize(num_poly_points);
                    }
            };

            struct ProcessState {
                size_t input_processed;
                size_t output_processed;
                std::vector<T> output;

                ProcessState() : input_processed(0), output_processed(0) {}
            };

            Denominator calculateDenominator(double rate, unsigned int num_poly_points, FirFilter<T, float>* filter = nullptr);
            bool canProcess(DenominatorImmutable* denom, DenominatorState* denomState, size_t readerAvailable, size_t writterAvailable, double rate);

            ProcessState process(DenominatorImmutable* denom, DenominatorState* denomState, const std::vector<T>& input, size_t readerAvailable, size_t writerWriteable, double rate);

        private:

    };

    template <typename T>
    class StereoFractionalDecimator: public Module<T, T> {
        
        public:
            enum class StereoDemodMode {
                SdrPll,
                AnalyticQuadrature
            };

            StereoFractionalDecimator(float rateMPX, float rate, float tau, unsigned int num_poly_points, FirFilter<T, float>* filter = nullptr);
            ~StereoFractionalDecimator();
            
            bool canProcess() override;
            void process() override;

            double deemphasisFilter(double sample, bool is_left_channel) {
                if (is_left_channel) {
                    deemph_state_L = (1.0 - deemph_alpha) * sample + deemph_alpha * deemph_state_L;
                    return deemph_state_L;
                } else {
                    deemph_state_R = (1.0 - deemph_alpha) * sample + deemph_alpha * deemph_state_R;
                    return deemph_state_R;
                }
            }

            // Calibrates the L-R subcarrier recovery gain at runtime. The default 2.0
            // assumes a standard FM-stereo MPX where the L-R subcarrier and (L+R) mono
            // signal share the same modulation index; if your transmitter uses different
            // coefficients you can fine-tune this to maximize stereo separation.
            void setStereoFactor(double factor) { stereo_factor = factor; }
            double getStereoFactor() const { return stereo_factor; }

            // Phase offset (in degrees) applied to the recovered 38 kHz subcarrier reference.
            //   0   -> cos(2ωt) — matches transmitters that put the L-R subcarrier in phase
            //                    with the squared 19 kHz pilot.
            //   90  -> sin(2ωt) — matches the standard FCC/EBU FM-stereo convention.
            // Sweep this to maximise channel separation if you don't know your TX's convention.
            void setSubcarrierPhase(double degrees) {
                subcarrier_phase_rad = degrees * (M_PI / 180.0);
                subcarrier_phase_cos = std::cos(subcarrier_phase_rad);
                subcarrier_phase_sin = std::sin(subcarrier_phase_rad);
            }
            double getSubcarrierPhase() const { return subcarrier_phase_rad * (180.0 / M_PI); }

            // Pilot envelope smoothing time constant (in seconds). Smaller values track
            // faster amplitude variations (good for unstable / fading transmitters);
            // larger values reject more pilot noise. Default = 3 ms, comfortably faster
            // than any realistic multipath fade rate.
            void setPilotEnvelopeTimeConstant(double seconds) {
                if (seconds <= 0.0) return;
                env_sq_alpha = 1.0 - std::exp(-1.0 / (inputSampleRate * seconds));
            }

            // Stereo↔mono blend thresholds, expressed in pilot RMS amplitude
            // (= √env_sq_smoothed, same scale as the MPX input).
            //   pilot_low_rms  : below this the decoder collapses to mono (no stereo noise)
            //   pilot_high_rms : above this full stereo is restored
            //   in between     : linear crossfade
            // Default initial values are 0.005 .. 0.020 (soft mono fallback enabled).
            // You can still disable manually by setting both to 0.
            void setStereoBlendThresholds(double pilot_low_rms, double pilot_high_rms) {
                blend_low_threshold  = pilot_low_rms;
                blend_high_threshold = pilot_high_rms;
            }
            double getStereoBlendLowThreshold()  const { return blend_low_threshold; }
            double getStereoBlendHighThreshold() const { return blend_high_threshold; }

            // Quiet-program de-hiss blend (audio-level dependent).
            // In very quiet passages we partially fold stereo toward mono to suppress
            // residual L-R noise/hiss. Full stereo is preserved once program level rises.
            void setQuietBlendThresholds(double audio_low, double audio_high, double min_stereo) {
                quiet_blend_low = audio_low;
                quiet_blend_high = audio_high;
                quiet_blend_min_stereo = min_stereo;
            }

            void setDemodMode(StereoDemodMode mode) { demod_mode = mode; }
            StereoDemodMode getDemodMode() const { return demod_mode; }

            // Fine tune MPX path delay (in input samples) to minimize stereo bleed.
            // In PLL mode this is often the most sensitive alignment knob.
            void setMpxAlignDelaySamples(size_t samples) {
                if (samples >= MPX_ALIGN_MAX_DELAY) samples = MPX_ALIGN_MAX_DELAY - 1;
                mpx_align_delay_samples = samples;
            }
            size_t getMpxAlignDelaySamples() const { return mpx_align_delay_samples; }
            
            private:
            unsigned int num_poly_points; //number of samples that the Lagrange interpolator will use
            float rate;
            FirFilter<T, float>* filter;
            
            // FMDemodMPX START
            unsigned int inputSampleRate;
            unsigned int outputSampleRate;

            void initializeFilters();
            bool initializedFilters_ = false;

            BiquadFilter* filter_19k;          // 19kHz pilot tone bandpass (subcarrier reference)
            MultistageFilter* filter_lp_lr;    // 15kHz LP, 8th order — applied to L-R baseband
            MultistageFilter* filter_lp_mono;  // 15kHz LP, 8th order — applied to mono baseband
                                               // (Both share identical coefficients; matched
                                               //  group delay and magnitude across the audio band.)
            // Pre-allocated work buffers (grown on demand, reused across process() calls)
            std::vector<double> input_fm_buf;
            std::vector<T> input_left_buf;
            std::vector<T> input_right_buf;

            // Runtime-tunable gain to convert recovered L-R baseband back to its original amplitude.
            // 2.0 is correct for any MPX where (L+R) mono and (L-R) subcarrier share the same
            // modulation index (the standard case).
            double stereo_factor;
            // Slow automatic trim on top of stereo_factor (user knob). Effective L-R gain is
            // stereo_factor * stereo_factor_auto.
            double stereo_factor_auto;
            double stereo_factor_auto_alpha;

            // Instantaneous envelope tracker for the analytic pilot pair (pilot, pilot_q).
            // env_sq_smoothed → A(t)² with a fast (~few-ms) IIR. Used to normalise the
            // recovered cos(2ωt) and sin(2ωt) references to unit amplitude on a
            // sample-by-sample basis. Replaces the older slow DC-tracker, which
            // struggled with amplitude-varying pilots (multipath fading, weak/unstable
            // transmitters).
            double env_sq_smoothed;
            double env_sq_alpha;

            // Optional stereo↔mono blend driven by the smoothed pilot RMS amplitude.
            //   pilot RMS ≤ blend_low_threshold  → pure mono (lr is muted)
            //   pilot RMS ≥ blend_high_threshold → full stereo
            //   in between                       → linear crossfade
            // Defaults: 0.005 .. 0.020 → blending enabled for weak/unstable pilots.
            double blend_low_threshold;
            double blend_high_threshold;
            double pilot_blend_smoothed;
            double pilot_blend_alpha;

            // Additional "quiet passage" stereo blend. Uses a smoothed |mono| envelope:
            // below quiet_blend_low -> keep only quiet_blend_min_stereo of L-R,
            // above quiet_blend_high -> full L-R, in between linear crossfade.
            double quiet_audio_env;
            double quiet_audio_alpha;
            double quiet_blend_low;
            double quiet_blend_high;
            double quiet_blend_min_stereo;

            // Quadrature-pilot generation via fractional delay (~2.526 samples @ 192 kHz / 19 kHz).
            // Lets us synthesise sin(2ωt) in addition to cos(2ωt), and therefore hit any
            // arbitrary subcarrier phase by mixing them.
            static const size_t PILOT_HISTORY_LEN = 8;
            double pilot_history[PILOT_HISTORY_LEN];
            size_t pilot_history_idx;
            size_t quad_delay_int;
            double quad_delay_frac;

            // Small MPX alignment delay (SDR++-style path matching).
            // Both mono and L-R branches consume the same delayed MPX sample, which
            // improves timing consistency against the pilot-derived reference path.
            static const size_t MPX_ALIGN_MAX_DELAY = 32;
            double mpx_delay_line[MPX_ALIGN_MAX_DELAY];
            size_t mpx_delay_index;
            size_t mpx_align_delay_samples;

            // Cached cos/sin of the user-selected subcarrier phase offset.
            double subcarrier_phase_rad;
            double subcarrier_phase_cos;
            double subcarrier_phase_sin;

            // PLL-based pilot reference generator (used for coherent 38 kHz NCO).
            PilotPLL* pilot_pll;
            // Smoothed polarity tracker between PLL and analytic 38 kHz references.
            // Keeps PLL path from occasionally settling with inverted stereo polarity.
            double pll_polarity_metric;
            double pll_polarity_alpha;
            int pll_polarity_sign;
            bool pll_polarity_locked;
            size_t pll_polarity_lock_samples;

            // Pilot RMS diagnostic metric (from the recovered pilot envelope).
            double pilot_strength;

            // Deemphasis filter states (50 µs time constant)
            double deemph_tau;
            double deemph_alpha;
            double deemph_state_L;
            double deemph_state_R;

            // Simple de-esser state: split deemphasized audio into low/high bands
            // with a 1-pole LP, then softly compress only the high band.
            double deesser_lp_a;
            double deesser_lp_L;
            double deesser_lp_R;

            // Slow per-channel DC blocker (~3 Hz @ inputSampleRate)
            double left_dc_offset, right_dc_offset;
            double balance_alpha;

            StereoDemodMode demod_mode;

            // FMDemodMPX STOP

            typename MonoFractionalDecimator<T>::DenominatorImmutable* denomImmutable;

            typename MonoFractionalDecimator<T>::DenominatorState* denomState_left;
            typename MonoFractionalDecimator<T>::DenominatorState* denomState_right;

            MonoFractionalDecimator<T> left_decimator;
            MonoFractionalDecimator<T> right_decimator;
    };

}