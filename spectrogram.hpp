#ifndef INCLUDE_SPECTROGRAM_HPP
#define INCLUDE_SPECTROGRAM_HPP

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
        Grid(float w, float h, int ws, int hs, bool is_log) : width(w), height(h),
            width_of_strips(w / ws), width_segments(ws), height_segments(hs), is_logarithmic(is_log);

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
            Strip(float w, float h, int s, bool is_log) : width(w), height(h), segments(s),
                segment_height(h / s), is_logarithmic(is_log);

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
        bool is_logarithmic
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
        //  sample_buffer_size: determines how many past samples the Spectrogram should store
        //  fft_window_size: determines how many bins the Short-Time Fourier Transform uses
        //  w,h: the width and height of the spectrogram respectively
        //  is_log: switch for linear or logarithmic frequency display
        Spectrogram(int samplerate, int fft_window_size, float x, float y, float w, float h, float ler, bool is_log) :
            bufferSize(samplerate), numBins(fft_window_size), grid_width(w* ler), legend_ratio(ler),
            stft(fft_window_size, fft_window_size / 4, 0, gam::HANN, gam::MAG_FREQ),
            grid(w* ler, h, samplerate, fft_window_size, is_log), width(w), height(h), x_offset(x), y_offset(y), is_logarithmic(is_log);

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

        int bufferSize, numBins;
        float width, height, x_offset, y_offset;
        float grid_width, legend_ratio;
        bool is_logarithmic;

        int bufferWrite;
        vector<vector<float>> spectrum;
        gam::STFT stft;
        Grid::Strip* legend_strip;
        float ref;
    };
}

#endif