
#include <cstdio>  // for printing to stdout
#include <algorithm>
#include <functional>

#include "spectrogram.hpp"

// using namespace gam;
using namespace std;

namespace al {
    Grid::Strip::Strip(float w, float h, int s, bool is_log) : width(w), height(h), segments(s), segment_height(h / s), is_logarithmic(is_log)
    {
        zero_color = getColor(0.f);
        float segment_scale_factor = static_cast<float>(segments) / log10(segments);
        this->primitive(TRIANGLE_STRIP);
        for (int i = 0; i <= segments; i++) {
            this->color(zero_color);
            this->color(zero_color);
            if (is_logarithmic) {
                float log_normalized = log10(static_cast<float>(i) + 1.f) / log10(segments + 1.f);
                float y_pos_log = height * log_normalized - (height / 2);
                this->vertex(-width / 2, y_pos_log);
                this->vertex(width / 2, y_pos_log); 
            }
            else {
                float y_pos_lin = segment_height * static_cast<float>(i) - height / 2;
                this->vertex(-width / 2, y_pos_lin);
                this->vertex(width / 2, y_pos_lin);
            }
        }
    }

    void Grid::Strip::map_value_to_color(vector<float> vals) {
        // need to pass in vector of values between zero and one
        // zero is for no amplitude, one is for max
        // first value is the lowest pitch bucket, last is highest
        // vector needs to have size 'segments' or 'segments + 1'
        for (int i = 0; i < segments; i++) {
            RGB curCol = get_color(vals[i]);
            this->colors()[2 * i] = curCol;
            this->colors()[2 * i + 1] = curCol;
        }
    }

    RGB Grid::Strip::get_color(float val) {
        // ripped from matlab's jet color map
        // 0.1242 0.3747 0.6253 0.8758
        if (val < 0) val = 0;
        if (val > 1) val = 1;
        float dr, dg, db;
        if (val < 0.1242) {
            db = 0.504 + ((1. - 0.504) / 0.1242) * val;
            dr = dg = 0;
        }
        else if (val < 0.3747) {
            db = 1.;
            dr = 0.;
            dg = (val - 0.1242) * (1. / (0.3747 - 0.1242));
        }
        else if (val < 0.6253) {
            db = (0.6253 - val) * (1. / (0.6253 - 0.3747));
            dg = 1.;
            dr = (val - 0.3747) * (1. / (0.6253 - 0.3747));
        }
        else if (val < 0.8758) {
            db = 0.;
            dr = 1.;
            dg = (0.8758 - val) * (1. / (0.8758 - 0.6253));
        }
        else {
            db = 0.;
            dg = 0.;
            dr = 1. - (val - 0.8758) * ((1. - 0.504) / (1. - 0.8758));
        }
        return RGB(dr, dg, db);
    }

    
    Grid::Grid(float w, float h, int ws, int hs, bool is_log) : 
        width(w), height(h), width_of_strips(w / ws),
        width_segments(ws), height_segments(hs), is_logarithmic(is_log) 
    {
        int strip_scale_factor = 1;
        for (int i = 0; i < width_segments; i++) {
            Strip* thisStrip = new Strip(width_of_strips, height, height_segments / strip_scale_factor, is_logarithmic);
            strips.push_back(thisStrip);
        }
        writePointer = 0;
    }

    Grid::~Grid() {
        for (int i = 0; i < width_segments; i++) {
            delete strips[i];
        }
    }

    void Grid::write_data(vector<float> vals) {
        strips[writePointer]->map_value_to_color(vals);
        writePointer = (writePointer + 1) % width_segments;
    }

    Grid::Strip* Grid::read_data(int index) {
        return strips[(writePointer + index) % width_segments];
    }
    
    Spectrogram::Spectrogram(int samplerate, int fft_window_size, float x, float y, float w, float h, float ler, bool is_log) :
        bufferSize(samplerate), numBins(fft_window_size), grid_width(w* ler), legend_ratio(ler),
        stft(fft_window_size, fft_window_size / 4, 0, gam::HANN, gam::MAG_FREQ),
        grid(w* ler, h, samplerate, fft_window_size, is_log), width(w), height(h), x_offset(x), y_offset(y), is_logarithmic(is_log)
    {
        spectrum.resize(bufferSize);
        for (int i = 0; i < bufferSize; i++) {
            spectrum[i].resize(numBins);
        }
        bufferWrite = 0;
        ref = 0.50;
        InitLegend();
    }
    
    void Spectrogram::write_sample(float sample) {
        if (stft(sample)) {
            for (int i = 0; i < numBins; ++i) {
                float val = stft.bin(i).real();
                spectrum[bufferWrite][i] = 1.5 * (log10(32 * (val)+ref) + 0.3);
            }
            grid.write_data(spectrum[bufferWrite]);
            bufferWrite = (bufferWrite + 1) % bufferSize;
        }
    }
    
    void Spectrogram::draw(Graphics& g) {
        // draw spectrogram
        g.loadIdentity();
        g.meshColor();
        g.translate(x_offset, y_offset);
        g.translate((-width + grid_width / bufferSize) / 2.f, height / 2);
        for (int i = 0; i < bufferSize; i++) {
            g.draw(*grid.read_data(i));
            g.translate(grid_width / bufferSize, 0);
        }
        // draw legend
        g.loadIdentity();
        g.translate(x_offset, y_offset);
        g.translate((width - (1 - legend_ratio) * 1.5 * width / 2) / 2.f, height / 2);
        g.draw(*legend_strip);
    }
    
    void Spectrogram::process_audio(AudioIOData& io) {
        while (io()) {
            write_sample(io.out(0));
        }
    }

    void Spectrogram::InitLegend() {
        vector<float> gradient;
        for (int i = 0; i < numBins; i++) {
            gradient.push_back(i * (height / 2) / numBins);
        }
        legend_strip = new Grid::Strip((1 - legend_ratio) * 1.5 * width / 2, height, numBins, false);
        legend_strip->map_value_to_color(gradient);
    }
}