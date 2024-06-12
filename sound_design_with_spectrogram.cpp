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

#include "src/spectrogram.hpp"

using namespace std;
using namespace al;

class AcidTest : public SynthVoice{
public:
    //gam::Accum<> tmr;
    //gam::Saw<> saw;
    gam::Square<> square;
    gam::Biquad<> lpf;
    gam::ADSR<> lpf_env;
    gam::ADSR<> amp_env;
    gam::Reson<> a; 
    //gam::OnePole<> freq;
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

        //tmr.freq(120. / 60. * 4.);
        //tmr.phaseMax();

        //freq.lag(0.1);

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
    Spectrogram mySpec{ 200, 4096, 512, 0, 0, 1.9, 1.9, 0.9, true};

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
        imguiBeginFrame();
        synthManager.drawSynthControlPanel();
        imguiEndFrame();
    } 

    void onDraw(Graphics &g) override { // Draw function
        g.clear();
        g.camera(Viewpoint::IDENTITY);
        mySpec.draw(g);
        imguiDraw();
    }

    void onSound(AudioIOData &io) override { // Audio callback
        synthManager.render(io);
        //mySpec.process_audio(io);

        // note that this 'while' can only happen in one place
        while (io()) { 
            mySpec.write_sample(io.out(0));
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
