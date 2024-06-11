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
    class Grid : public Mesh{
    public:
        // Constructor for class Grid
        //  w : width (0 - 2) of grid
        //  h : height (0 - 2) of grid
        //  ws : number of horizontal segments - how many strips wide
        //  hs : number of vertical segments - how many segments high
        //  is_log : boolean switch for whether to distribute vertical segments linearly or logarithmically
        Grid(float r, float h, int ws, int hs, bool is_log) :
            radius(r), height(h), width_of_strips(r* M_2PI / ws),
            width_segments(ws), height_segments(hs), is_logarithmic(is_log)
        {
            this->primitive(TRIANGLES);
            RGB zero_color = get_color(0.f);
            this->vertices_bank = new Vertex * [width_segments + 1];
            this->colors_bank = new RGB * [width_segments + 1];

            float theta = 0;
            for (int i = 0; i <= width_segments; i++) {
                this->vertices_bank[i] = new Vertex[height_segments + 1];
                this->colors_bank[i] = new RGB[height_segments + 1];

                if (!this->colors_bank[i]) {
                    std::cerr << "Memory allocation for colors_bank[" << i << "] failed." << std::endl;
                    // Handle the error, possibly with cleanup and exit
                }
                float cartX = radius * cos(theta);
                float cartZ = radius * sin(theta);
                for (int j = 0; j <= height_segments; j++) {
                    this->colors_bank[i][j] = zero_color;
                    if (is_log) {
                        float log_normalized = log10((float)j + 1.f) / log10(height_segments + 1.f);
                        float z_pos_log = height * log_normalized - (height / 2);
                        this->vertices_bank[i][j] = Vertex(cartX, z_pos_log, cartZ);
                    }
                    else {
                        float z_pos_lin = (height_segments / height) * (float)j - (height / 2);
                        this->vertices_bank[i][j] = Vertex(cartX, z_pos_lin, cartZ);
                    }
                }
                theta -= width_of_strips / radius;
            }

            for (int i = 0; i < width_segments; i++) {
                for (int j = 0; j < height_segments; j++) {
                    this->vertex(this->vertices_bank[i][j]);
                    this->color(this->colors_bank[i][j]);
                    this->vertex(this->vertices_bank[i + 1][j]);
                    this->color(this->colors_bank[i][j]);
                    this->vertex(this->vertices_bank[i + 1][j + 1]);
                    this->color(this->colors_bank[i][j]);
                    this->vertex(this->vertices_bank[i][j]);
                    this->color(this->colors_bank[i][j]);
                    this->vertex(this->vertices_bank[i][j + 1]);
                    this->color(this->colors_bank[i][j]);
                    this->vertex(this->vertices_bank[i + 1][j + 1]);
                    this->color(this->colors_bank[i][j]);
                }
            }
            //int strip_scale_factor = 1;
            //float theta = 0;
            //for (int i = 0; i < width_segments; i++) {
            //    theta -= width_of_strips / radius;
            //    Strip* thisStrip = new Strip(width_of_strips, height, radius, theta, height_segments / strip_scale_factor, is_logarithmic);
            //    strips.push_back(thisStrip);
            //}
            writePointer = 0;
        }

        /*
        * A strip is a mesh of customizable dimension that is two vertices wide and user specified vertices tall.
        * Strips are meant to display data using color, and currently use the same color scheme as Matlab's "Jet" color palette.
        */
        

        // Grid destructor
        ~Grid() {
            for (int i = 0; i <= width_segments; i++) {
                free(this->vertices_bank[i]);
                free(this->colors_bank[i]);
            }
            free(this->vertices_bank);
            free(this->colors_bank);
            //for (int i = 0; i < width_segments; i++) {
            //    delete strips[i];
            //}
        }

        // writes a std::vector of floats 
        void write_data(const std::vector<float>& vals) {
            for (int j = 0; j < height_segments; j++) {
                RGB curCol = get_color(vals[j]);
                int offset = ((writePointer * height_segments) + j) * 6;
                for (int k = 0; k < 6; k++) {
                    this->colors()[offset + k] = curCol;
                }

            }
            writePointer = (writePointer + 1) % width_segments;
        }

    protected:
        RGB get_color(float val) {
            // ripped from matlab's jet color map
            // 0.1242 0.3747 0.6253 0.8758
            if (val < 0) val = 0;
            if (val > 1) val = 1;
            float dr, dg, db;
            if (val < 0.1242) {
                db = val * 3.9936 + 0.5042;
                //db = 0.504 + ((1. - 0.504) / 0.1242) * val;
                dr = dg = 0;
            }
            else if (val < 0.3747) {
                db = 1.;
                dr = 0.;
                dg = val * 3.9920 - 0.4958;
            }
            else if (val < 0.6253) {
                db = 2.495 - val * 3.9904;
                dg = 1.;
                dr = val * 3.990 - 1.4952;
            }
            else if (val < 0.8758) {
                db = 0.;
                dr = 1.;
                dg = 3.4962 - val * 3.9920;
            }
            else {
                db = 0.;
                dg = 0.;
                dr = 4.4962 - val * 3.9920;
            }
            return RGB(dr, dg, db);
        }
        Vertex** vertices_bank;
        RGB** colors_bank;
        float radius, height, width_of_strips;
        int width_segments, height_segments;
        int writePointer;
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
        /**
         * @brief 
         * 
         
         */
        struct SpectrogramSphericalParams{
            int buffer_size;
            int fft_window_size;
            int hop_size;
            float w;
            float h;
            bool is_log;
            bool is_sphere;
            // is_sphere params
            float radius;
        };

        Spectrogram(SpectrogramSphericalParams);

        // Constructor for class Spectrogram
        //  buffer_size: determines how many past STFT outputs the Spectrogram should store
        //  fft_window_size: determines how many bins the Short-Time Fourier Transform uses
        //  x,z: the x and y coordinate offset for the spectrogram display
        //  w,h: the width and height of the spectrogram respectively
        //  ler: legend ratio - what percentage of the width and height for spectrogram vs legend
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
        //void InitLegend();

        int bufferSize, numBins, hopSize;
        float width, height, x_offset, z_offset;
        float grid_width, legend_ratio;
        float radius;
        bool is_logarithmic, is_spherical;

        float spherical_angle;

        int bufferWrite;
        std::vector<std::vector<float>> spectrum;
        gam::STFT stft;
        float ref;
    };



    // Spectrogram flat constructor
    Spectrogram::Spectrogram(int samplerate, int fft_window_size, int hop_size, float x, float z, float w, float h, float ler, bool is_log) :
        bufferSize(samplerate), numBins(fft_window_size), hopSize(hop_size), grid_width(w* ler), legend_ratio(ler),
        stft(fft_window_size, hop_size, 0, gam::HANN, gam::MAG_FREQ),
        grid(w* ler, h, samplerate, fft_window_size, is_log), width(w), height(h), x_offset(x), z_offset(z), spherical_angle(0.f),
        is_logarithmic(is_log)//, is_spherical(is_sphere)
    {
        spectrum.resize(bufferSize);
        for (int i = 0; i < bufferSize; i++) {
            spectrum[i].resize(numBins);
        }
        bufferWrite = 0;
        ref = 0.50;
        //InitLegend();
    }

    // Spectrogram spherical constructor
    Spectrogram::Spectrogram(SpectrogramSphericalParams params) :
        bufferSize(params.buffer_size), numBins(params.fft_window_size), hopSize(params.hop_size), grid_width(params.w * M_PI), legend_ratio(0.0f),
        stft(params.fft_window_size, params.hop_size, 0, gam::HANN, gam::MAG_FREQ),
        grid(params.w * M_PI, params.h, params.buffer_size, params.fft_window_size, params.is_log), width(params.w * M_PI), height(params.h), x_offset(0), z_offset(0),
        is_logarithmic(params.is_log), is_spherical(params.is_sphere), radius(params.radius), spherical_angle((params.w * M_PI / params.buffer_size) / params.radius)
    {
        spectrum.resize(bufferSize);
        for (int i = 0; i < bufferSize; i++) {
            spectrum[i].resize(numBins);
        }
        bufferWrite = 0;
        ref = 0.50;
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
        // g.loadIdentity();
        g.meshColor();
        if(is_spherical){
            //g.draw(*grid.read_data(0));
            g.draw(grid);
            /*for (int i = 0; i < bufferSize; i++) {
                g.draw(*grid.read_data(i));
            }*/
        } 
        // else {
        //     g.translate(x_offset, z_offset);
        //     //if(is_spherical) g.translate()
        //     g.translate((-width + grid_width / bufferSize) / 2.f, 0);
        //     for (int i = 0; i < bufferSize; i++) {
        //         g.draw(*grid.read_data(i));
        //         g.translate(grid_width / bufferSize, 0);
        //     }
        //     // draw legend
        //     g.loadIdentity();
        //     g.translate(x_offset, z_offset);
        //     g.translate((width - (1 - legend_ratio) * 1.5 * width / 2) / 2.f, 0);
        //     g.draw(*legend_strip);
        // }
    }

    void Spectrogram::process_audio(AudioIOData& io) {
        while (io()) {
            write_sample(io.out(0));
        }
    }

    /*void Spectrogram::InitLegend() {
        std::vector<float> gradient;
        for (int i = 0; i < numBins; i++) {
            gradient.push_back(i * (height / 2) / numBins);
        }
        legend_strip = new Grid::Strip((1 - legend_ratio) * 1.5 * width / 2, height, numBins, false);
        legend_strip->map_value_to_color(gradient);
    }*/
}

#endif