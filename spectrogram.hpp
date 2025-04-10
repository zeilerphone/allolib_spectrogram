#ifndef INCLUDE_SPECTROGRAM_HPP
#define INCLUDE_SPECTROGRAM_HPP

#include <algorithm>
#include <functional>

#include "Gamma/Analysis.h"
#include "Gamma/Gamma.h"
#include "Gamma/DFT.h"

#include "al/graphics/al_Shapes.hpp"
#include "al/types/al_Color.hpp"

namespace al {
    /*
    A Grid manages writing to and reading from a vector of Strip objects.
    The number of strip objects and parameters for the strip objects are parameters in the Grid constructor.
    */    
    class Grid {
    public:
        // Constructor for class Grid
        //  w : width (0 - 2) of grid
        //  h : height (0 - 2) of grid
        //  ws : number of horizontal segments - how many strips wide
        //  hs : number of vertical segments - how many segments high
        //  is_log : boolean switch for whether to distribute vertical segments linearly or logarithmically
        Grid(float w, float h, int ws, int hs, bool is_log);

        /*
        * A strip is a mesh of customizable dimension that is two vertices wide and user specified vertices tall.
        * Strips are meant to display data using color, and currently use the same color scheme as Matlab's "Jet" color palette.
        */
        class Strip : public Mesh {
        public:
            // Constructor for class Strip
            //  w : width (2 is screen width) of strip
            //  h : height (2 is screen height) of strip
            //  s : number of segments the strip should be split into
            //  is_log : boolean switch for whether to distribute strips linearly or logarithmically
            Strip(float w, float h, int s, bool is_log);

            /* 
                Assigns strip vertex colors based off of input std::vector<float>.
                The size of the vector must be equal to or less than `segments`. 
                The values in the vector must lie within the range [0,1].
            */
            void map_value_to_color(std::vector<float> vals);

        protected:
            RGB get_color(float val);
            float width, height, segment_height;
            int segments;
            RGB zero_color;
            bool is_logarithmic;
        };

        // Grid destructor
        ~Grid();

        // writes a std::vector of floats into the next Strip
        void write_data(std::vector<float> vals);

        // returns a pointer to the Strip that is 'index' places past the current Strip
        Strip* read_data(int index);

    protected:
        float width, height, width_of_strips;
        int width_segments, height_segments;
        int writePointer;
        std::vector<Strip*> strips;
        bool is_logarithmic;
    };

    /*
    A Spectogram displays the frequency decomposition of a digital signal over time. The x-axis is time,
    the y-axis is frequency, and the color represents the amplitude of each frequency.
    The Spectrogram class consists of a Grid object, a legend to define which colors correspond to which freqencies,
    a Short-Time Fourier transform object from Gamma, and internal variables and parameters.
    */
    class Spectrogram {
    public:
        // Constructor for class Spectrogram
        //  buffer_size: determines how many past STFT outputs the Spectrogram should store and display
        //  fft_window_size: determines how many bins the Short-Time Fourier Transform uses
        //  hop_size: determines how many samples the STFT should wait for before calculating the next set of spectrum values
        //  x,z: the x and y coordinate offset for the spectrogram display
        //  w,h: the width and height of the spectrogram respectively
        //  ler: legend ratio - what percentage of the width for spectrogram vs legend
        //  is_log: switch for linear or logarithmic frequency display
        Spectrogram(int buffer_size, int fft_window_size, int hop_size, float x, float z, float w, float h, float ler, bool is_log);

        // writes a sample from audio output to the internal STFT and buffers
        void write_sample(float sample);

        // draws the Spectrogram on graphics object 'g'
        void draw(Graphics& g);

        // records audio data from channel 0 to the internal STFT and buffers
        // only use and pass 'io' if this is the only thing using the audio data
        // otherwise use Spectrogram::writeSample instead
        void process_audio(AudioIOData& io);

        Grid grid;

    protected:
        void InitLegend();

        int bufferSize, numBins, hopSize;
        float width, height, x_offset, z_offset;
        float grid_width, legend_ratio;
        bool is_logarithmic;

        int bufferWrite;
        std::vector<std::vector<float>> spectrum;
        gam::STFT stft;
        Grid::Strip* legend_strip;
        float ref;
    };



    Grid::Strip::Strip(float w, float h, int s, bool is_log) : width(w), height(h), segments(s), segment_height(h / s), is_logarithmic(is_log)
    {
        zero_color = get_color(0.f);
        float segment_scale_factor = static_cast<float>(segments) / log10(segments);
        this->primitive(TRIANGLE_STRIP);
        for (int i = 0; i <= segments; i++) {
            this->color(zero_color);
            this->color(zero_color);
            if (is_logarithmic) {
                float log_normalized = log10(static_cast<float>(i) + 1.f) / log10(segments + 1.f);
                float z_pos_log = height * log_normalized - (height / 2);
                this->vertex(-width / 2, z_pos_log);
                this->vertex(width / 2, z_pos_log);
            }
            else {
                float z_pos_lin = segment_height * static_cast<float>(i) - height / 2;
                this->vertex(-width / 2, z_pos_lin);
                this->vertex(width / 2, z_pos_lin);
            }
        }
    }

    void Grid::Strip::map_value_to_color(std::vector<float> vals) {
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

    void Grid::write_data(std::vector<float> vals) {
        strips[writePointer]->map_value_to_color(vals);
        writePointer = (writePointer + 1) % width_segments;
    }

    Grid::Strip* Grid::read_data(int index) {
        return strips[(writePointer + index) % width_segments];
    }

    Spectrogram::Spectrogram(int samplerate, int fft_window_size, int hop_size, float x, float z, float w, float h, float ler, bool is_log) :
        bufferSize(samplerate), numBins(fft_window_size), hopSize(hop_size), grid_width(w* ler), legend_ratio(ler),
        stft(fft_window_size, hop_size, 0, gam::HANN, gam::MAG_FREQ),
        grid(w* ler, h, samplerate, fft_window_size, is_log), width(w), height(h), x_offset(x), z_offset(z), 
        is_logarithmic(is_log)
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
        g.translate(x_offset, z_offset);
        g.translate((-width + grid_width / bufferSize) / 2.f, 0);
        for (int i = 0; i < bufferSize; i++) {
            g.draw(*grid.read_data(i));
            g.translate(grid_width / bufferSize, 0);
        }
        // draw legend
        g.loadIdentity();
        g.translate(x_offset, z_offset);
        g.translate((width - (1 - legend_ratio) * 1.5 * width / 2) / 2.f, 0);
        g.draw(*legend_strip);
        g.loadIdentity();
    }

    void Spectrogram::process_audio(AudioIOData& io) {
        while (io()) {
            write_sample(io.out(0));
        }
    }

    void Spectrogram::InitLegend() {
        std::vector<float> gradient;
        for (int i = 0; i < numBins; i++) {
            gradient.push_back(i * (height / 2) / numBins);
        }
        legend_strip = new Grid::Strip((1 - legend_ratio) * 1.5 * width / 2, height, numBins, false);
        legend_strip->map_value_to_color(gradient);
    }
}

#endif