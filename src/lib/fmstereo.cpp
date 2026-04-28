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

#include "fmstereo.hpp"

#include <algorithm>
#include <vector>
#include <complex>
#include <cmath>
#include <iostream>

#define TEST_DIRECTFMINPUT false

using namespace Csdr;

template <typename T>
MonoFractionalDecimator<T>::MonoFractionalDecimator() {

}

template <typename T>
typename MonoFractionalDecimator<T>::Denominator MonoFractionalDecimator<T>::calculateDenominator(float rate, unsigned int num_poly_points, FirFilter<T, float> *filter) {

    Denominator denom;

    denom.num_poly_points = (num_poly_points &~ 1);
    denom.poly_precalc_denomiator.resize(denom.num_poly_points);
    denom.xifirst = (-(denom.num_poly_points / 2) + 1);
    denom.xilast = (denom.num_poly_points / 2);
    denom.coeffs_buf.resize(denom.num_poly_points);
    denom.rate = (rate);

    int id = 0; //index in poly_precalc_denomiator
    for (int xi = denom.xifirst; xi <= denom.xilast; xi++) {
        denom.poly_precalc_denomiator[id] = 1;
        for(int xj = denom.xifirst; xj <= denom.xilast; xj++) {
            //poly_precalc_denomiator could be integer as well. But that would later add a necessary conversion.
            if (xi != xj) denom.poly_precalc_denomiator[id] *= (xi - xj);
        }
        id++;
    }

    denom.where = -denom.xifirst;

    return denom;
}

template <typename T>
MonoFractionalDecimator<T>::~MonoFractionalDecimator() {
    
}

template <typename T>
bool MonoFractionalDecimator<T>::canProcess(DenominatorImmutable* denomImmutable, DenominatorState* denomState, size_t readerAvailable, size_t writterAvailable, float rate) {   
    size_t size = std::min(readerAvailable, (size_t) ceilf((writterAvailable) / rate));
    size_t filterLen = denomImmutable->filter != nullptr ? denomImmutable->filter->getOverhead() : 0;
    return ceilf(denomState->where) + denomImmutable->num_poly_points + filterLen < size;
}

template <typename T>
typename MonoFractionalDecimator<T>::ProcessState MonoFractionalDecimator<T>::process(DenominatorImmutable* denom, DenominatorState* denomState, std::vector<T> input, size_t readerAvailable, size_t writerWriteable, float rate) {
    
    if (!denom || !denomState) {
        printf("MonoFractionalDecimator::ProcessState::process - Denominators not available!\n");
        return ProcessState(); // Return empty state
    }
    
    ProcessState state;
    state.output = std::vector<T>(writerWriteable);

    int oi = 0; //output index
    int index_high, index;
    size_t size = std::min(readerAvailable, (size_t) ceilf(writerWriteable / rate));
    size_t filterLen = denom->filter != nullptr ? denom->filter->getOverhead() : 0;

    int ceil_val = static_cast<int>(ceilf(denomState->where));
    index_high = ceil_val;

    int num_poly = denom->num_poly_points;
    int limit = ceil_val + num_poly + filterLen;

    //we optimize to calculate ceilf(where) only once every iteration, so we do it here:
    while (limit < size) {
        // num_poly_points above is theoretically more than we could have here, but this makes the spectrum look good
        index = index_high - 1;

        int id = 0;
        float xwhere = denomState->where - index;
        for (int xi = denom->xifirst; xi <= denom->xilast; xi++) {
            denomState->coeffs_buf[id] = 1;
            for (int xj = denom->xifirst; xj <= denom->xilast; xj++) {
                if (xi != xj) denomState->coeffs_buf[id] *= (xwhere - xj);
            }
            id++;
        }
        
        T acc = T(0);
        if (denom->filter != nullptr) {
            for (int i = 0; i < denom->num_poly_points; i++) {
                SparseView<float> sparse = denom->filter->sparse(input.data());
                acc += (denomState->coeffs_buf[i] / denom->poly_precalc_denomiator[i]) * sparse[index + i];
            }
        } else {
            for (int i = 0; i < denom->num_poly_points; i++) {
                acc += (denomState->coeffs_buf[i] / denom->poly_precalc_denomiator[i]) * (input[index + i]);
            }
        }

        state.output[oi++] = acc;
        denomState->where += rate;

        ceil_val = static_cast<int>(ceilf(denomState->where));
        index_high = ceil_val;

        num_poly = denom->num_poly_points;
        limit = ceil_val + num_poly + filterLen;
    }

    int input_processed = 0;

    input_processed = index + denom->xifirst;
    denomState->where -= input_processed;

    state.input_processed = input_processed;
    state.output_processed = oi;

    return state;
}

template <typename T>
void StereoFractionalDecimator<T>::initializeFilters() {
    std::lock_guard<std::mutex> lock(this->processMutex);

    filter_19k     = new BiquadFilter();
    filter_lp_lr   = new MultistageFilter();
    filter_lp_mono = new MultistageFilter();
    // Keep PLL loop bandwidth moderate so pilot lock remains robust.
    pilot_pll      = new PilotPLL(inputSampleRate, 19000.0, 0.707, 20.0);

    // 19 kHz pilot bandpass — feeds the 38 kHz reference generator.
    // Wider BW reduces time-domain ringing on this implementation and avoids
    // HF artifacts observed with overly narrow pilot filtering.
    filter_19k->setBandpass2(19000.0, 1800.0, inputSampleRate);
    // Identical 14.2 kHz LP @ 8th order on both mono and L-R paths so that magnitude
    // and group delay match exactly across the audio band — the key for clean
    // stereo separation. The 16th-order Butterworth gives:
    //    ~−33 dB at 19 kHz  (pilot residue → effectively inaudible)
    //    ~−59 dB at 23 kHz  (lower L-R subcarrier sideband — main crosstalk source at 8th order)
    //    ~−143 dB at 38 kHz (subcarrier centre)
    // The extra 4 biquads per path are cheap on modern CPUs and give a roughly 30 dB
    // improvement in stop-band rejection over the 8th-order baseline.
    // Keep full FM audio brightness around the traditional 15 kHz limit.
    filter_lp_lr  ->setLowpass(15000.0, inputSampleRate, 8);
    filter_lp_mono->setLowpass(15000.0, inputSampleRate, 8);
    // TODO: make it adjustable
    // Deemphasis time constant (50 microseconds)
    deemph_alpha = exp(-(1.0 / inputSampleRate) / deemph_tau);
    deemph_state_L = 0.0;
    deemph_state_R = 0.0;
    // 1-pole LP split around ~3.5 kHz for de-essing.
    deesser_lp_a = std::exp(-2.0 * M_PI * 3500.0 / inputSampleRate);
    deesser_lp_L = 0.0;
    deesser_lp_R = 0.0;

    // Slow per-channel DC blocker (~3 Hz) — protects against tiny carrier bias.
    left_dc_offset = right_dc_offset = 0.0;
    balance_alpha  = 0.0001;

    // Conservative default for wider transmitter compatibility. Some custom TX chains
    // over-drive the decoded L-R path with factor 2.0, which can manifest as HF artifacts.
    // You can still tune at runtime via setStereoFactor().
    stereo_factor = 1.5;

    // Instantaneous envelope normaliser: we use p² + p_q² ≈ A(t)² to scale the
    // recovered references on every sample, so that varying pilot amplitudes don't
    // modulate the audio. ~3 ms IIR smoothing knocks down the residual 38 kHz
    // ripple from the imperfect quadrature pilot while still tracking realistic
    // (multipath / weak-signal) fade rates of 100 Hz or so.
    env_sq_smoothed = 0.0;
    env_sq_alpha    = 1.0 - std::exp(-1.0 / (inputSampleRate * 0.003));

    // Gentle stereo fallback for real-air reception: when pilot gets weak/noisy,
    // fade L-R toward mono instead of letting high-band artifacts break through.
    blend_low_threshold  = 0.004;
    blend_high_threshold = 0.014;
    pilot_blend_smoothed = 0.0;
    pilot_blend_alpha    = 1.0 - std::exp(-1.0 / (inputSampleRate * 0.08)); // ~80 ms

    // Quiet-program de-hiss blend disabled by default for maximum linearity.
    // Re-enable via setQuietBlendThresholds() if you want low-level hiss masking.
    quiet_audio_env = 0.0;
    quiet_audio_alpha = 1.0 - std::exp(-1.0 / (inputSampleRate * 0.05)); // ~50 ms
    quiet_blend_low = 0.0;
    quiet_blend_high = 0.0;
    quiet_blend_min_stereo = 1.0;

    // Quadrature pilot via fractional delay: π/2 phase shift @ 19 kHz corresponds to
    // exactly fs / (4·19000) samples. We split into integer + fractional parts and
    // interpolate linearly between two history samples.
    double quad_delay_samples = inputSampleRate / (4.0 * 19000.0);
    quad_delay_int  = static_cast<size_t>(std::floor(quad_delay_samples));
    quad_delay_frac = quad_delay_samples - quad_delay_int;
    pilot_history_idx = 0;
    for (size_t i = 0; i < PILOT_HISTORY_LEN; ++i) pilot_history[i] = 0.0;

    // Default to sin-form reference (90° offset) — this matches the standard FCC/EBU
    // FM-stereo MPX where the L-R subcarrier is (L−R)·sin(2ωt). For non-standard
    // transmitters that use the cos convention call setSubcarrierPhase(0.0).
    subcarrier_phase_rad = M_PI / 2.0;
    subcarrier_phase_cos = 0.0;
    subcarrier_phase_sin = 1.0;

    initializedFilters_ = true;
}

template <typename T>
StereoFractionalDecimator<T>::StereoFractionalDecimator(float rateMPX, float rate, float tau, unsigned int num_poly_points, FirFilter<T, float> *filter):
    num_poly_points(num_poly_points &~ 1),
    inputSampleRate(rateMPX), outputSampleRate(rateMPX),
    rate(rate),
    deemph_tau(tau),
    filter(nullptr)
{
    try {
        initializeFilters();

        pilot_strength = 0.0;

        denomImmutable = new typename MonoFractionalDecimator<T>::DenominatorImmutable(this->num_poly_points, rate, filter);
        denomState_left = new typename MonoFractionalDecimator<T>::DenominatorState(-denomImmutable->xifirst, denomImmutable->num_poly_points);
        denomState_right = new typename MonoFractionalDecimator<T>::DenominatorState(-denomImmutable->xifirst, denomImmutable->num_poly_points);
        
    }
    catch (const std::exception& e) {
        printf("StereoFractionalDecimator::Constructor exception: %s\n", e.what());
        throw;
    }
}

template <typename T>
StereoFractionalDecimator<T>::~StereoFractionalDecimator() {
    delete denomState_right;
    delete denomState_left;
    delete denomImmutable;

    delete filter_19k;
    delete filter_lp_lr;
    delete filter_lp_mono;
    delete pilot_pll;
    deemph_state_L = deemph_state_R = 0.0;
    pilot_strength = 0.0;
    left_dc_offset = right_dc_offset = 0.0;
}

template <typename T>
bool StereoFractionalDecimator<T>::canProcess() {
    std::lock_guard<std::mutex> lock(this->processMutex);
    
    size_t size = std::min(this->reader->available(), (size_t) ceilf((this->writer->writeable() / 2) / rate));
    size_t filterLen = denomImmutable->filter != nullptr ? denomImmutable->filter->getOverhead() : 0;
    return ceilf(denomState_left->where) + denomImmutable->num_poly_points + filterLen < size;
}

template <typename T>
void StereoFractionalDecimator<T>::process() {
    std::lock_guard<std::mutex> lock(this->processMutex);
    
    int oi = 0;
    int input_frames_processed = 0;
    
    size_t available_samples = this->reader->available();
    size_t available_stereo_frames = available_samples / 2;
    size_t writeable_samples = this->writer->writeable();
    size_t writeable_stereo_frames = writeable_samples / 2;
    size_t size_frames = std::min(available_stereo_frames, (size_t) std::ceil(writeable_stereo_frames / rate));
    size_t filterLen = filter != nullptr ? filter->getOverhead() : 0;
 
    T* input = this->reader->getReadPointer();
    T* output = this->writer->getWritePointer();

    typename MonoFractionalDecimator<T>::ProcessState state_left;
    typename MonoFractionalDecimator<T>::ProcessState state_right;

    bool leftCanProcess = left_decimator.canProcess(denomImmutable, denomState_left, available_stereo_frames, writeable_stereo_frames, rate);

    if(leftCanProcess)
    {
        if(!initializedFilters_){
            printf("StereoFractionalDecimator::process - Filters not initialized!\n");
            std::terminate();
        }

        // Grow work buffers on demand only; reuse across calls to avoid per-call heap churn.
        if (input_left_buf.size() < available_samples) {
            input_fm_buf.resize(available_samples);
            input_left_buf.resize(available_samples);
            input_right_buf.resize(available_samples);
        }

        for(size_t i = 0; i < available_samples; i++) {
            input_fm_buf[i] = input[i];
        }

        size_t i_left = 0;
        size_t i_right = 0;

        // Pure FM-stereo decoder loop:
        //   mono = LP(mpx)
        //   lr   = LP(mpx * coherent_38k * stereo_factor)
        //
        // The 38 kHz reference is generated from the filtered pilot and its
        // quadrature counterpart, with per-sample envelope normalisation.
        // This keeps channel assignment stable while remaining robust against
        // pilot-amplitude fluctuations.
        //
        // Both audio paths run through identical 16th-order Butterworth low-passes
        // (15 kHz), so their amplitude and group delay match exactly across the
        // audio band — the key for clean, frequency-independent separation.
        for(size_t i = 0; i < available_samples; i++) {
            double mpx_signal = input_fm_buf[i];

            // 19 kHz pilot extraction.
            double pilot_19k = filter_19k->process(mpx_signal);

            // Slip-free 38 kHz reference from pilot + quadrature pilot.
            // Unlike PLL lock, this path cannot jump by 180° in the audio matrix,
            // so left/right channels won't intermittently swap.
            pilot_history[pilot_history_idx] = pilot_19k;
            size_t i_newer = (pilot_history_idx + PILOT_HISTORY_LEN - quad_delay_int    ) % PILOT_HISTORY_LEN;
            size_t i_older = (pilot_history_idx + PILOT_HISTORY_LEN - quad_delay_int - 1) % PILOT_HISTORY_LEN;
            double pilot_q = (1.0 - quad_delay_frac) * pilot_history[i_newer]
                           + quad_delay_frac        * pilot_history[i_older];
            pilot_history_idx = (pilot_history_idx + 1) % PILOT_HISTORY_LEN;

            double env_sq_inst = pilot_19k * pilot_19k + pilot_q * pilot_q;
            env_sq_smoothed   += env_sq_alpha * (env_sq_inst - env_sq_smoothed);

            double cos_38k = 0.0;
            double sin_38k = 0.0;
            if (env_sq_smoothed > 1e-12) {
                cos_38k = (pilot_19k * pilot_19k - pilot_q * pilot_q) / env_sq_smoothed;
                sin_38k = -2.0 * (pilot_19k * pilot_q) / env_sq_smoothed;
            }

            // Tunable phase: cos(2ωt − θ) = cos(θ)·cos(2ωt) + sin(θ)·sin(2ωt).
            // Calibrate via setSubcarrierPhase() for your transmitter convention.
            double coherent_38k = subcarrier_phase_cos * cos_38k
                                + subcarrier_phase_sin * sin_38k;

            pilot_strength = std::sqrt(env_sq_smoothed);

            // Mono and L-R baseband, both filtered by identical 15 kHz LPs.
            double mono = filter_lp_mono->process(mpx_signal);
            double lr   = filter_lp_lr  ->process(mpx_signal * coherent_38k * stereo_factor);

            // Optional stereo↔mono blend: gradually mute lr below the user-defined
            // pilot-strength window. Disabled by default (low == high == 0).
            if (blend_high_threshold > blend_low_threshold) {
                pilot_blend_smoothed += pilot_blend_alpha * (pilot_strength - pilot_blend_smoothed);
                double blend = (pilot_blend_smoothed - blend_low_threshold)
                             / (blend_high_threshold - blend_low_threshold);
                if (blend < 0.0) blend = 0.0;
                if (blend > 1.0) blend = 1.0;
                lr *= blend;
            }

            // Quiet passage de-hiss: reduce L-R noise floor when program level is very low.
            quiet_audio_env += quiet_audio_alpha * (std::fabs(mono) - quiet_audio_env);
            if (quiet_blend_high > quiet_blend_low) {
                double q = (quiet_audio_env - quiet_blend_low)
                         / (quiet_blend_high - quiet_blend_low);
                if (q < 0.0) q = 0.0;
                if (q > 1.0) q = 1.0;
                double quiet_blend = quiet_blend_min_stereo + (1.0 - quiet_blend_min_stereo) * q;
                lr *= quiet_blend;
            }

            double left_raw  = mono + lr;
            double right_raw = mono - lr;

            // Slow DC blocker (~3 Hz @ inputSampleRate) — only removes numerical/carrier
            // bias, doesn't touch audio.
            left_dc_offset  += balance_alpha * (left_raw  - left_dc_offset);
            right_dc_offset += balance_alpha * (right_raw - right_dc_offset);
            left_raw  -= left_dc_offset;
            right_raw -= right_dc_offset;

            // Apply deemphasis first (broadcast FM standard), then a very gentle
            // soft-limiter with extra headroom. This significantly reduces
            // sibilant/speech harshness versus clipping before deemphasis.
            double left_deemph  = deemphasisFilter(left_raw,  true);
            double right_deemph = deemphasisFilter(right_raw, false);

            // De-esser: split into low/high with a 1-pole LP, then softly compress
            // only the high band. This keeps clarity while reducing "sz/s" harshness.
            deesser_lp_L = deesser_lp_a * deesser_lp_L + (1.0 - deesser_lp_a) * left_deemph;
            deesser_lp_R = deesser_lp_a * deesser_lp_R + (1.0 - deesser_lp_a) * right_deemph;
            double high_L = left_deemph  - deesser_lp_L;
            double high_R = right_deemph - deesser_lp_R;

            auto deess_high = [](double h) -> double {
                const double threshold = 0.030;
                const double ratio = 0.35; // compress above threshold
                double ah = std::fabs(h);
                if (ah <= threshold) return h;
                double sign = (h >= 0.0) ? 1.0 : -1.0;
                double compressed = threshold + (ah - threshold) * ratio;
                return sign * compressed;
            };

            double left_out  = deesser_lp_L + deess_high(high_L);
            double right_out = deesser_lp_R + deess_high(high_R);

            // Global headroom + transparent peak guard to avoid overall "overdriven"
            // character after stereo matrixing/de-essing.
            const double output_gain = 0.64;
            auto peak_guard = [](double x) -> double {
                if (x > 0.98) return 0.98;
                if (x < -0.98) return -0.98;
                return x;
            };
            input_left_buf [i_left++]  = peak_guard(output_gain * left_out);
            input_right_buf[i_right++] = peak_guard(output_gain * right_out);
        }

#if !TEST_DIRECTFMINPUT
        state_left  = left_decimator.process(denomImmutable,  denomState_left,  input_left_buf,  i_left,  writeable_stereo_frames, rate);
        state_right = right_decimator.process(denomImmutable, denomState_right, input_right_buf, i_right, writeable_stereo_frames, rate);

        size_t max_input_consumed = (state_left.input_processed + state_right.input_processed) / 2;
        size_t min_output = state_left.output_processed + state_right.output_processed;

        for (size_t i = 0; i < min_output / 2; i++) {
            output[i * 2]     = state_left.output[i];
            output[i * 2 + 1] = state_right.output[i];
        }

        size_t input_samples_consumed  = max_input_consumed;
        size_t output_samples_produced = min_output;
#else
        size_t max_input_consumed = (i_left + i_right) / 2;
        size_t min_output = i_left + i_right;

        for (size_t i = 0; i < i_left; i++) {
            output[i * 2]     = input_left_buf[i];
            output[i * 2 + 1] = input_right_buf[i];
        }

        size_t input_samples_consumed  = max_input_consumed;
        size_t output_samples_produced = min_output;
#endif

        this->reader->advance(input_samples_consumed);
        this->writer->advance(output_samples_produced);
    }
}

namespace Csdr {
    template class MonoFractionalDecimator<float>;
    template class StereoFractionalDecimator<float>;
}