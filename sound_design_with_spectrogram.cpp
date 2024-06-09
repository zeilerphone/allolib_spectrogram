#include <cstdio>  // for printing to stdout
#include <algorithm>
#include <functional>

#include "al/app/al_App.hpp"

#include "Gamma/Oscillator.h"
#include "Gamma/Effects.h"
#include "Gamma/Filter.h"
#include "Gamma/Envelope.h"
#include "Gamma/Types.h"
#include "Gamma/DFT.h"

#include "al/scene/al_PolySynth.hpp"
#include "al/scene/al_SynthSequencer.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"

using namespace std;
using namespace al;

//class System1 : public SynthVoice {
//
//};

//class Oscilliscope : public Mesh {
//public:
//    Oscilliscope(int samplerate) : bufferSize(samplerate) {
//        this->primitive(Mesh::LINE_STRIP);
//        for (int i = 0; i < bufferSize; i++) {
//            this->vertex((i / static_cast<float>(bufferSize)) * 1.8f - 0.9f, 0);
//            buffer.push_back(0.f);
//        }
//        bufferWrite = 0;
//    }
//
//    void writeSample(float sample) {
//        buffer[bufferWrite] = sample;
//        bufferWrite = (bufferWrite + 1) % bufferSize;
//    }
//
//    void update() {
//        for (int i = 0; i < bufferSize; i++) {
//            this->vertices()[i][1] = buffer[(bufferWrite + i) % bufferSize];
//        }
//    }
//
//protected:
//    int bufferSize;
//    int bufferWrite;
//    vector<float> buffer;
//};

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
    Strip(float w, float h, int s, bool is_log) : width(w), height(h), segments(s), segment_height(h / s), is_logarithmic(is_log)
    {
        zero_color = getColor(0.);
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

    // assigns strip vertex colors based off of input std::vector 
    // size must match
    void map_value_to_color(vector<float> vals) {
        // need to pass in vector of values between zero and one
        // zero is for no amplitude, one is for max
        // first value is the lowest pitch bucket, last is highest
        // vector needs to have size 'segments' or 'segments + 1'
        for (int i = 0; i < segments; i++) {
            RGB curCol = getColor(vals[i]);
            this->colors()[2 * i] = curCol;
            this->colors()[2 * i + 1] = curCol;
        }
    }

protected:
    RGB getColor(float val) {
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
    float width, height, segment_height;
    int segments;
    RGB zero_color;
    bool is_logarithmic;
};

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
    Grid(float w, float h, int ws, int hs, bool is_log) :
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

    // Grid destructor
    ~Grid() {
        for (int i = 0; i < width_segments; i++) {
            delete strips[i];
        }
    }

    // writes a std::vector of floats into the next Strip
    void write_data(vector<float> vals) {
        strips[writePointer]->map_value_to_color(vals);
        writePointer = (writePointer + 1) % width_segments;
    }

    // returns a pointer to the Strip that is 'index' places past the current Strip
    Strip* read_data(int index) {
        return strips[(writePointer + index) % width_segments];
    }

protected:
    float width, height, width_of_strips;
    int width_segments, height_segments;
    int writePointer;
    vector<Strip*> strips;
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
    //  sample_buffer_size: determines how many past samples the Spectrogram should store
    //  fft_window_size: determines how many bins the Short-Time Fourier Transform uses
    //  x,y: the x and y coordinate offset for the spectrogram display
    //  w,h: the width and height of the spectrogram respectively
    //  ler: legend ratio - what percentage of the width and height for spectrogram vs legend
    //  is_log: switch for linear or logarithmic frequency display
    Spectrogram(int samplerate, int fft_window_size, float x, float y, float w, float h, float ler, bool is_log) :
        bufferSize(samplerate), numBins(fft_window_size), grid_width(w * ler), legend_ratio(ler), 
        stft(fft_window_size, fft_window_size / 4, 0, gam::HANN, gam::MAG_FREQ),
        grid(w * ler, h, samplerate, fft_window_size, is_log), width(w), height(h), x_offset(x), y_offset(y), is_logarithmic(is_log)
    {
        spectrum.resize(bufferSize);
        for (int i = 0; i < bufferSize; i++) {
            spectrum[i].resize(numBins);
        }
        bufferWrite = 0;
        ref = 0.50;

        InitLegend();
    }

    // writes a sample from audio output to the internal STFT and buffers
    void writeSample(float sample) {
        if (stft(sample)) {
            for (int i = 0; i < numBins; ++i) {
                float val = stft.bin(i).real();
                spectrum[bufferWrite][i] = 1.5 * (log10(32 * (val)+ref) + 0.3);
            }
            grid.write_data(spectrum[bufferWrite]);
            bufferWrite = (bufferWrite + 1) % bufferSize;
        }
    }

    // draws the Spectrogram on graphics object 'g'
    void draw(Graphics& g) {
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
        // draw labels (unimplemented)
    }

    // records audio data from channel 0 to the internal STFT and buffers
    // only use and pass 'io' if this is the only thing using the audio data
    // otherwise use Spectrogram::writeSample instead
    void process_audio(AudioIOData& io) {
        while (io()) {
            writeSample(io.out(0));
        }
    }
    Grid grid;

protected:
    void InitLegend() {
        vector<float> gradient;
        for (int i = 0; i < numBins; i++) {
            gradient.push_back(i * (height / 2) / numBins);
        }
        legend_strip = new Strip((1 - legend_ratio) * 1.5 * width / 2, height, numBins, false);
        legend_strip->map_value_to_color(gradient);
    }
    int bufferSize, numBins;
    float width, height, x_offset, y_offset;
    float grid_width, legend_ratio;
    bool is_logarithmic;

    int bufferWrite;
    vector<vector<float>> spectrum;
    gam::STFT stft;
    Strip* legend_strip;
    float ref;
};

class AcidTest : public SynthVoice{
public:
    gam::Accum<> tmr;
    gam::Saw<> saw;
    gam::Square<> square;
    gam::Biquad<> lpf;
    gam::ADSR<> lpf_env;
    gam::ADSR<> amp_env;
    gam::OnePole<> freq;
    gam::Pan<> pan;
    float ref, clip;
    //int seqLength;

    float sample, last_sample;
    float lpf_env_depth;

    void init() override {
        lpf.type(gam::LOW_PASS);

        lpf.res(20);
        lpf.freq(440);

        amp_env.attack(0.01);
        amp_env.decay(0.5);
        amp_env.sustain(0.5);
        amp_env.release(0.01);

        lpf_env.attack(0.01);
        lpf_env.decay(1);
        lpf_env.sustain(0.5);
        lpf_env.release(0.5);
        //lpf_env.amp(0.75f);

        tmr.freq(120. / 60. * 4.);
        tmr.phaseMax();

        freq.lag(0.1);

        ref = 27.5; // A0 = 27.5 Hz
        //int seqLength = 8;
        //pitches = new float[8] { 11, 11, 23, 11, 11, 21, 6, 11 };
        //int seqLength = 16;
        lpf_env_depth = 0.75;

        createInternalTriggerParameter("Amp", 0.3, 0.0, 1.0);
        createInternalTriggerParameter("f", 60, 20, 5000);
        createInternalTriggerParameter("keyboard_octave", 2, 0, 8);

        createInternalTriggerParameter("LPF:resonance", 4.5, 0.01, 10.0);
        createInternalTriggerParameter("LPF:frequency_offset", 2, 0.01, 5);

        createInternalTriggerParameter("AmpEnv:Attack",    0.05,   0.01,   3.0);
        createInternalTriggerParameter("AmpEnv:Decay",     0.4,    0.01,   3.0);
        createInternalTriggerParameter("AmpEnv:Sustain",   0.5,    0.01,   3.0);
        createInternalTriggerParameter("AmpEnv:Release",   0.2,    0.1,    5.0);

        createInternalTriggerParameter("LPFFreqEnv:Attack",    0.05,   0.01,   3.0);
        createInternalTriggerParameter("LPFFreqEnv:Decay",     0.25,   0.01,   3.0);
        createInternalTriggerParameter("LPFFreqEnv:Sustain",   0.01,   0.01,   3.0);
        createInternalTriggerParameter("LPFFreqEnv:Release",   0.2,    0.1,    5.0);
        createInternalTriggerParameter("LPFFreqEnv:Depth",     64,     0.01,   128);

        createInternalTriggerParameter("pan", 0.0, -1.0, 1.0);

        clip = 0.5f;
    }

    virtual void onProcess(AudioIOData& io) override {
        updateFromParameters();
        float amp = getInternalParameterValue("Amp");
        float LPFFreqOffset = getInternalParameterValue("LPF:frequency_offset");
        float frequency = getInternalParameterValue("f");
        while (io()) {
            sample = square();                          // get the next sample from the square wave 
            float amp_e = amp_env();                    // get the next value from the amplitude 
            float lpf_e = lpf_env() * lpf_env_depth;    // get the next value from the low-pass filter envelope and multiply by the LPF depth parameter
            float lpf_frequency = (LPFFreqOffset + lpf_e) * frequency;  // calculate the frequency the low pass filter 
            lpf.freq(lpf_frequency);    // set the LPF's frequency to calculated frequency (this changing quickly along with resonance causes clicking noise)

            /*
                The LPF's frequency response involves unity gain (no amplitude change) for frequencies below the cutoff, and zero gain for frequencys past it.
                The resonance of the LPF causes the amplitude of frequencies at the cutoff frequency to be amplified above unity gain (positive amplitude change).
                Because the LPF's frequency is modulated, and the default attack for the frequency modulation envelope is low, the LPF's frequency sweeps quickly and causes a lot of different frequencies to be amplified at roughly the same time, which causes a "clicking" noise

                In the future, this can be fixed by reducing the resonance to zero during the attack of the LPF's envelope
            */

            sample = lpf(sample) * amp_e * amp;         // pass the sample through the lowpass filter and multiply it by the amplitude envelope and amplitude parameter
            
            clip = 2 * amp;     // set the clip to be twice the amplitude
            
            // clip the sample based on it
            if (sample >  clip) sample =  clip;
            if (sample < -clip) sample = -clip;

            float s1;
            float s2;
            
            pan(sample, s1, s2);    // pass the sample through the pan object
            io.out(0) += s1;        // send the samples to the respective audio output channels
            io.out(1) += s2;
            last_sample = sample;
        }
        // when the amplitude envelope is finished, stop generating the signal
        if (amp_env.done()) free();
    }

    virtual void onTriggerOn() override {
        updateFromParameters();
        amp_env.reset();
        lpf_env.reset();
    }

    virtual void onTriggerOff() override {
        amp_env.triggerRelease();
        lpf_env.triggerRelease();
    }

    void updateFromParameters() {
        float frequency = getInternalParameterValue("f");
        square.freq(frequency);

        lpf.res(getInternalParameterValue("LPF:resonance"));
        lpf.freq(getInternalParameterValue("LPF:frequency_offset") * frequency);

        amp_env.attack(getInternalParameterValue("AmpEnv:Attack"));
        amp_env.decay(getInternalParameterValue("AmpEnv:Decay"));
        amp_env.sustain(getInternalParameterValue("AmpEnv:Sustain"));
        amp_env.release(getInternalParameterValue("AmpEnv:Release"));

        lpf_env.attack(getInternalParameterValue("LPFFreqEnv:Attack"));
        lpf_env.decay(getInternalParameterValue("LPFFreqEnv:Decay"));
        lpf_env.sustain(getInternalParameterValue("LPFFreqEnv:Sustain"));
        lpf_env.release(getInternalParameterValue("LPFFreqEnv:Release"));
        lpf_env_depth = getInternalParameterValue("LPFFreqEnv:Depth");
    }
};

class MyApp: public App {
public:
    SynthGUIManager<AcidTest> synthManager{ "" };
    Spectrogram mySpec{ 800, 512, 0, 0, 1.9, 1.9, 0.9, false };
    //Oscilliscope myScope{ 44100 };

    void onInit() override { // Called on app start
        std::cout << "onInit()" << std::endl;
        imguiInit();
        navControl().active(false);
    }

    void onCreate() override { // Called when graphics context is available
        std::cout << "onCreate()" << std::endl;
        synthManager.synthSequencer().playSequence("AcidTest.synthSequence");
    }

    void onAnimate(double dt) override { // Called once before drawing
        /*color += 0.01;
        if (color > 1.0) {
            color -= 1.0;
        }*/
        //oscillator.freq(220 + 880 * color);
        imguiBeginFrame();
        synthManager.drawSynthControlPanel();
        imguiEndFrame();
    } 

    void onDraw(Graphics &g) override { // Draw function
        //myScope.update();
        
        g.clear();
        g.camera(Viewpoint::IDENTITY);
        mySpec.draw(g);
        g.loadIdentity();

        //g.color(1);
        //g.draw(myScope);
        imguiDraw();
    }

    void onSound(AudioIOData &io) override { // Audio callback
        synthManager.render(io);
        //mySpec.process_audio(io);

        // note that this 'while' can only happen in one place
        while (io()) { 
            mySpec.writeSample(io.out(0));
            //myScope.writeSample(io.out(0));
        }
    }

    void onMessage(osc::Message &m) override { // OSC message callback
        m.print();
    }

    bool onKeyDown(Keyboard const& k) override {
        int midiNote = asciiToMIDI(k.key());
        if (midiNote > 0) {
            float oct = floor(synthManager.voice()->getInternalParameter("keyboard_octave"));
            synthManager.voice()->setInternalParameterValue(
                "f", ::pow(2.f, (midiNote - 69.f) / 12.f) * 27.5f * oct);
            synthManager.triggerOn(midiNote);
        }
        return true;
    }

    bool onKeyUp(Keyboard const& k) override {
        int midiNote = asciiToMIDI(k.key());
        if (midiNote > 0) {
            synthManager.triggerOff(midiNote);
        }
        return true;
    }

    void onExit() override { imguiShutdown(); }
};


int main()
{
    MyApp app;
    auto dev = AudioDevice("devicename");
    app.configureAudio(dev, 44100, 256, 2, 2);

    app.start();
    return 0;
}
