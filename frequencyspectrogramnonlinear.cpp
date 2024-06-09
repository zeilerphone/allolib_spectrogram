
#include <cstdio>  // for printing to stdout
#include <algorithm>
#include <functional>

#include "Gamma/Analysis.h"
#include "Gamma/Effects.h"
#include "Gamma/Envelope.h"
#include "Gamma/Gamma.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Types.h"
#include "Gamma/DFT.h"

#include "al/app/al_App.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/scene/al_PolySynth.hpp"
#include "al/scene/al_SynthSequencer.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"

// using namespace gam;
using namespace al;
using namespace std;

// tables for oscillator
gam::ArrayPow2<float>
tbSaw(2048), tbSqr(2048), tbImp(2048), tbSin(2048), tbPls(2048),
tb__1(2048), tb__2(2048), tb__3(2048), tb__4(2048);

//template <typename T>
//std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
//{
//    assert(a.size() == b.size());
//
//    std::vector<T> result;
//    result.reserve(a.size());
//
//    std::transform(a.begin(), a.end(), b.begin(),
//        std::back_inserter(result), std::plus<T>());
//    return result;
//}
//
//template <typename T>
//std::vector<T> operator/(const std::vector<T>& a, const T b)
//{
//    std::vector<T> result;
//    std::vector<T> bvect(a.size(), b);
//    result.reserve(a.size());
//
//    std::transform(a.begin(), a.end(), bvect.begin(),
//        std::back_inserter(result), std::divides<T>());
//    return result;
//}

// Oscilliscope that inherits from mesh 
//class Oscilliscope : public Mesh {
//public:
//    Oscilliscope(int samplerate) : bufferSize(samplerate) {
//        this->primitive(Mesh::LINE_STRIP);
//        for (int i = 0; i < bufferSize; i++) {
//            this->vertex((i / static_cast<float>(bufferSize)) * 2.f - 1.f, 0);
//            buffer.push_back(0.f);
//        }
//    }
//
//    void writeSample(float sample) {
//        for (int i = 0; i < bufferSize - 1; i++) {
//            buffer[i] = buffer[i + 1];
//        }
//        buffer[bufferSize - 1] = sample;
//    }
//
//    void update() {
//        for (int i = 0; i < bufferSize; i++) {
//            this->vertices()[i][1] = buffer[i];
//        }
//    }
//
//protected:
//    int bufferSize;
//    vector<float> buffer;
//};

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
        for (int i = 0; i < strips.size(); i++) {
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

class Spectrogram {
public:
    // Constructor for class Spectrogram
    //  sample_buffer_size: determines how many past samples the Spectrogram should store
    //  fft_window_size: determines how many bins the Short-Time Fourier Transform uses
    //  w,h: the width and height of the spectrogram respectively
    //  ler: legend ratio - what percentage of the width and height for spectrogram vs legend
    //  lar: label ratio - unimplemented
    //  is_log: switch for linear or logarithmic frequency display
    Spectrogram(int samplerate, int fft_window_size, float w, float h, float ler, float lar, bool is_log) :
        bufferSize(samplerate), numBins(fft_window_size), grid_width(w* ler), legend_ratio(ler), label_ratio(lar),
        stft(fft_window_size, fft_window_size / 4, 0, gam::HANN, gam::MAG_FREQ),
        grid(w* ler, h* lar, samplerate, fft_window_size, is_log), width(w), height(h), is_logarithmic(is_log)
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

            //grid.write_data(spectrum[bufferWrite]);
            bufferWrite = (bufferWrite + 1) % bufferSize;
        }
        /*spectrum_LTR[LTRWrite] = spectrum_LTR[LTRWrite] + (spectrum[bufferWrite] / static_cast<float>(time_resolution_factor));
        time_resolution_counter++;
        if (time_resolution_counter == time_resolution_factor) {
            time_resolution_counter = 0;
            grid.write_data(spectrum_LTR[LTRWrite]);
            LTRWrite = (LTRWrite + 1) % bufferSize / time_resolution_factor;
        }*/

    }

    // draws the Spectrogram on graphics object 'g'
    void draw(Graphics& g) {
        // draw spectrogram
        g.loadIdentity();
        g.meshColor();
        g.translate((-width + grid_width / bufferSize) / 2.f, -(1 - label_ratio) * height / 2);
        for (int i = 0; i < bufferSize; i++) {
            g.draw(*grid.read_data(i));
            g.translate(grid_width / bufferSize, 0);
        }
        // draw legend
        g.loadIdentity();
        g.translate((width - (1 - legend_ratio) * 1.5 * width / 2) / 2.f, -(1 - label_ratio) * height / 2);
        g.draw(*legend_strip);
        // draw labels
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
        legend_strip = new Strip((1 - legend_ratio) * 1.5 * width / 2, label_ratio * height, numBins, false);
        legend_strip->map_value_to_color(gradient);
    }

    int bufferSize;
    int numBins;
    int bufferWrite;
    vector<vector<float>> spectrum;
    gam::STFT stft;
    Strip* legend_strip;
    float ref;
    float width, height;
    float grid_width;
    float legend_ratio, label_ratio;
    bool is_logarithmic;
};


class FrequencyScope : public Mesh {
public:
    FrequencyScope(int fft_window_size) :
        numBins(fft_window_size),
        stft(fft_window_size, fft_window_size / 4, 0, gam::HANN, gam::MAG_FREQ)
    {
        this->reset();
        this->primitive(Mesh::LINE_STRIP);
        for (int i = 0; i < numBins; i++) {
            this->vertex(i / static_cast<float>(numBins) * 1.9f - 0.95f, -0.3);
        }
        spectrum.resize(numBins);
        ref = 0.5;
    }

    void writeSample(float sample) {
        if (stft(sample)) {
            for (int i = 0; i < numBins; ++i) {
                // cout << stft.bin(i).real() << endl;
                //spectrum[i] = pow(stft.bin(i).real() * 10, 1.3) * 5;
                float val = stft.bin(i).real();
                //if (val > max) max = val;
                //spectrum[i] =  20 * log10(32.f * val - ref) - 0.f;
                spectrum[i] = log10(32 * (val) + ref) + 0.3;
                //cout << val << ", " << max << endl;
                // 0.03
            }
        }
    }

    void update() {
        //this->reset();
        //this->primitive(Mesh::LINE_STRIP);
        for (int i = 0; i < numBins; i++) {
            //this->vertex(i / static_cast<float>(numBins) * 4.f - 1.f, spectrum[i], 0.0);
            this->vertices()[i][1] = spectrum[i];
        }
    }
protected:
    int numBins;
    vector<float> spectrum;
    gam::STFT stft;
    float ref;
    float max = 0;
};

class OscTrm : public SynthVoice {
public:

    // Unit generators
    gam::Pan<> mPan;
    gam::Sine<> mTrm;
    gam::Osc<> mOsc;
    gam::ADSR<> mTrmEnv;
    //gam::Env<2> mTrmEnv;
    gam::ADSR<> mAmpEnv;
    gam::EnvFollow<> mEnvFollow;  // envelope follower to connect audio output to graphics

    // Additional members
    Mesh mMesh;

    // Initialize voice. This function will nly be called once per voice
    virtual void init() {

        // Intialize envelope
        mAmpEnv.curve(0); // make segments lines
        mAmpEnv.levels(0, 0.3, 0.3, 0); // These tables are not normalized, so scale to 0.3
        //        mAmpEnv.sustainPoint(1); // Make point 2 sustain until a release is issued

        mTrmEnv.curve(0);
        mTrmEnv.levels(0, 1, 1, 0);
        //        mTrmEnv.sustainPoint(1); // Make point 2 sustain until a release is issued

                // We have the mesh be a sphere
        addDisc(mMesh, 1.0, 30);

        createInternalTriggerParameter("amplitude", 0.3, 0.0, 1.0);
        createInternalTriggerParameter("frequency", 60, 20, 5000);
        createInternalTriggerParameter("attackTime", 0.1, 0.01, 3.0);
        createInternalTriggerParameter("releaseTime", 3.0, 0.1, 10.0);
        createInternalTriggerParameter("sustain", 0.7, 0.0, 1.0);
        createInternalTriggerParameter("curve", 4.0, -10.0, 10.0);
        createInternalTriggerParameter("pan", 0.0, -1.0, 1.0);
        createInternalTriggerParameter("table", 0, 0, 8);
        createInternalTriggerParameter("trm1", 3.5, 0.2, 20);
        createInternalTriggerParameter("trm2", 5.8, 0.2, 20);
        createInternalTriggerParameter("trmRise", 0.5, 0.1, 2);
        createInternalTriggerParameter("trmDepth", 0.1, 0.0, 1.0);
    }

    //
    virtual void onProcess(AudioIOData& io) override {
        //updateFromParameters();
        float oscFreq = getInternalParameterValue("frequency");
        float amp = getInternalParameterValue("amplitude");
        float trmDepth = getInternalParameterValue("trmDepth");
        while (io()) {

            mTrm.freq(mTrmEnv());
            //float trmAmp = mAmp - mTrm()*mTrmDepth; // Replaced with line below
            float trmAmp = (mTrm() * 0.5 + 0.5) * trmDepth + (1 - trmDepth); // Corrected
            float s1 = mOsc() * mAmpEnv() * trmAmp * amp;
            float s2;
            mEnvFollow(s1);
            mPan(s1, s1, s2);
            io.out(0) += s1;
            io.out(1) += s2;
        }
        // We need to let the synth know that this voice is done
        // by calling the free(). This takes the voice out of the
        // rendering chain
        if (mAmpEnv.done() && (mEnvFollow.value() < 0.001f)) free();
    }

    virtual void onProcess(Graphics& g) {
        float frequency = getInternalParameterValue("frequency");
        float amplitude = getInternalParameterValue("amplitude");
        g.pushMatrix();
        g.translate(amplitude, amplitude, -4);
        //g.scale(frequency/2000, frequency/4000, 1);
        float scaling = getInternalParameterValue("trmDepth");
        g.scale(scaling * frequency / 200, scaling * frequency / 400, scaling * 1);
        g.color(mEnvFollow.value(), frequency / 1000, mEnvFollow.value() * 10, 0.4);
        g.draw(mMesh);
        g.popMatrix();
    }

    virtual void onTriggerOn() override {
        updateFromParameters();

        mAmpEnv.reset();
        mTrmEnv.reset();

        // Map table number to table in memory
        switch (int(getInternalParameterValue("table"))) {
        case 0: mOsc.source(tbSaw); break;
        case 1: mOsc.source(tbSqr); break;
        case 2: mOsc.source(tbImp); break;
        case 3: mOsc.source(tbSin); break;
        case 4: mOsc.source(tbPls); break;
        case 5: mOsc.source(tb__1); break;
        case 6: mOsc.source(tb__2); break;
        case 7: mOsc.source(tb__3); break;
        case 8: mOsc.source(tb__4); break;
        }
    }

    virtual void onTriggerOff() override {
        mAmpEnv.triggerRelease();
        mTrmEnv.triggerRelease();
    }

    void updateFromParameters() {
        mOsc.freq(getInternalParameterValue("frequency"));
        mAmpEnv.attack(getInternalParameterValue("attackTime"));
        mAmpEnv.decay(getInternalParameterValue("attackTime"));
        mAmpEnv.release(getInternalParameterValue("releaseTime"));
        mAmpEnv.sustain(getInternalParameterValue("sustain"));
        mAmpEnv.curve(getInternalParameterValue("curve"));
        mPan.pos(getInternalParameterValue("pan"));

        mTrmEnv.levels(getInternalParameterValue("trm1"),
            getInternalParameterValue("trm2"),
            getInternalParameterValue("trm2"),
            getInternalParameterValue("trm1"));

        mTrmEnv.attack(getInternalParameterValue("trmRise"));
        mTrmEnv.decay(getInternalParameterValue("trmRise"));
        mTrmEnv.release(getInternalParameterValue("trmRise"));
    }
};


// We make an app.
class MyApp : public App
{
public:
    SynthGUIManager<OscTrm> synthManager{ "synth5" };
    //Oscilliscope myScope{ static_cast<int>(AudioIO().framesPerSecond()) };
    int scopeSampleRate = 800;
    int fft_window_size = 512;  // performs best with powers of 2 i.e. 256, 512, 1024
    bool is_logarithmic = true;
    Spectrogram mySpec{ 800, 512, 1.9, 1.9, 0.9, 1., true };
    //FrequencyScope myFreq{ 512 };
    //    ParameterMIDI parameterMIDI;
    int midiNote;


    virtual void onInit() override {
        imguiInit();
        navControl().active(false);  // Disable navigation via keyboard, since we
        // will be using keyboard for note triggering
// Set sampling rate for Gamma objects from app's audio
        gam::sampleRate(audioIO().framesPerSecond());


        gam::addSinesPow<1>(tbSaw, 9, 1);
        gam::addSinesPow<1>(tbSqr, 9, 2);
        gam::addSinesPow<0>(tbImp, 9, 1);
        gam::addSine(tbSin);

        {
            float A[] = { 1,1,1,1,0.7,0.5,0.3,0.1 };
            gam::addSines(tbPls, A, 8);
        }

        {
            float A[] = { 1, 0.4, 0.65, 0.3, 0.18, 0.08 };
            float C[] = { 1,4,7,11,15,18 };
            gam::addSines(tb__1, A, C, 6);
        }

        // inharmonic partials
        {
            float A[] = { 0.5,0.8,0.7,1,0.3,0.4,0.2,0.12 };
            float C[] = { 3,4,7,8,11,12,15,16 };
            gam::addSines(tb__2, A, C, 8);
        }

        // inharmonic partials
        {
            float A[] = { 1, 0.7, 0.45, 0.3, 0.15, 0.08 };
            float C[] = { 10, 27, 54, 81, 108, 135 };
            gam::addSines(tb__3, A, C, 6);
        }

        // harmonics 20-27
        {
            float A[] = { 0.2, 0.4, 0.6, 1, 0.7, 0.5, 0.3, 0.1 };
            gam::addSines(tb__4, A, 8, 20);
        }
    }
    void onCreate() override {
        // Play example sequence. Comment this line to start from scratch
        //    synthManager.synthSequencer().playSequence("synth2.synthSequence");
        synthManager.synthRecorder().verbose(true);
    }

    void onSound(AudioIOData& io) override {
        synthManager.render(io);  // Render audio
        mySpec.process_audio(io);
        while (io()) {
            //myScope.writeSample((io.out(0) + io.out(1)) / 2.f);
            //myFreq.writeSample(io.out(0));
        }
    }

    void onAnimate(double dt) override {
        imguiBeginFrame();
        synthManager.drawSynthControlPanel();
        imguiEndFrame();
        //myScope.update();
        //myFreq.update();
    }

    void onDraw(Graphics& g) override {
        //myFreq.update();
        g.clear();
        g.camera(Viewpoint::IDENTITY);
        mySpec.draw(g);
        //synthManager.render(g);
        g.loadIdentity();
        g.color(1);
        //g.draw(myFreq);
        // Draw GUI
        imguiDraw();
    }

    bool onKeyDown(Keyboard const& k) override {
        if (ParameterGUI::usingKeyboard()) {  // Ignore keys if GUI is using them
            return true;
        }
        if (k.shift()) {
            // If shift pressed then keyboard sets preset
            int presetNumber = asciiToIndex(k.key());
            synthManager.recallPreset(presetNumber);
        }
        else {
            // Otherwise trigger note for polyphonic synth
            int midiNote = asciiToMIDI(k.key());
            if (midiNote > 0) {
                synthManager.voice()->setInternalParameterValue(
                    "frequency", ::pow(2.f, (midiNote - 69.f) / 12.f) * 432.f);
                synthManager.triggerOn(midiNote);
            }
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

int main() {
    MyApp app;

    // Set up audio
    app.configureAudio(44100., 512, 2, 0);

    app.start();
}
