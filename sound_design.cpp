#include "al/app/al_App.hpp"
#include "Gamma/Oscillator.h"
#include "Gamma/Effects.h"
#include "Gamma/Filter.h"
#include "Gamma/Envelope.h"

#include "al/scene/al_PolySynth.hpp"
#include "al/scene/al_SynthSequencer.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"

using namespace al;

//class System1 : public SynthVoice {
//
//};

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
    int step;
    float* pitches;  
    //float* pitches2;  // D minor (ii)
    //float* pitches3;  // G major (V)
    //float* pitches4;  // C major (I)
    float ref;
    //int seqLength;
    float debugMax;

    float lpf_env_depth;

    void init() override {
        lpf.type(gam::LOW_PASS);

        lpf.res(20);
        lpf.freq(440);

        amp_env.attack(0.05);
        amp_env.decay(0.5);
        amp_env.sustain(0.5);
        amp_env.release(0.01);

        lpf_env.attack(0.05);
        lpf_env.decay(1);
        lpf_env.sustain(0.5);
        lpf_env.release(0.5);
        //lpf_env.amp(0.75f);

        tmr.freq(120. / 60. * 4.);
        tmr.phaseMax();

        freq.lag(0.1);

        step = 0;
        ref = 27.5; // A0 = 27.5 Hz
        //int seqLength = 8;
        //pitches = new float[8] { 11, 11, 23, 11, 11, 21, 6, 11 };
        //int seqLength = 16;
        pitches = new float[16] {13, 13, 21, 13, 13, 23, 13, 18, 13, 13, 21, 13, 25, 23, 18, 13  };
        lpf_env_depth = 0.75;

        createInternalTriggerParameter("Amp", 0.3, 0.0, 1.0);
        createInternalTriggerParameter("f", 60, 20, 5000);
        createInternalTriggerParameter("keyboard_octave", 2, 0, 8);

        createInternalTriggerParameter("LPF:resonance", 4.5, 0.01, 10.0);
        createInternalTriggerParameter("LPF:frequency_mult", 2, 0.01, 5);

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

        createInternalTriggerParameter("useSeq", 0, 0, 1);

        debugMax = 0;
    }

    virtual void onProcess(AudioIOData& io) override {
        updateFromParameters();
        float amp = getInternalParameterValue("Amp");
        float LPFFreqMult = getInternalParameterValue("LPF:frequency_mult");
        float useSequencer = getInternalParameterValue("useSeq");
        float frequency = getInternalParameterValue("f");
        while (io()) {
            if (abs(useSequencer - 1.f) <= 0.1) {
                if (tmr()) {
                    // ref was 50 -> about a G1, now A0
                    float f = ref * pow(2, (pitches[step]-1) / 12.);
                    //step = (step + 1) % 8;
                    step = (step + 1) % 16;
                    freq = f;
                   /* amp_env.resetSoft();
                    lpf_env.resetSoft();*/
                }
                square.freq(freq());
            }

            float s1 = square();

            float amp_e = amp_env();
            float lpf_e = lpf_env() * lpf_env_depth;
            float lpf_frequency = (LPFFreqMult + lpf_e) * frequency;
            
            lpf.freq(lpf_frequency);
            //lpf.freq(e * (modCutoff.paraU() * 6000 + 500) + 40);
            s1 = lpf(s1) * amp_e * amp;

            /*if (abs(s1) >= debugMax) {
                debugMax = abs(s1);
                std::cout << debugMax << std::endl;
            }*/

            float s2;
            pan(s1, s1, s2);
            io.out(0) += s1;
            io.out(1) += s2;
        }
        if (amp_env.done()) free();
    }

    virtual void onTriggerOn() override {
        updateFromParameters();
        amp_env.reset();
        lpf_env.reset();

    }

    virtual void onTriggerOff() override {
        amp_env.triggerRelease();
    }

    void updateFromParameters() {
        float frequency = getInternalParameterValue("f");
        square.freq(frequency);

        lpf.res(getInternalParameterValue("LPF:resonance"));
        lpf.freq(getInternalParameterValue("LPF:frequency_mult") * frequency);

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
    //AcidTest acid;
    //float color = 0.0;

   /* MyApp() : App(), acid() {
    }*/

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
        //g.clear(color);
        g.clear();
        synthManager.render(g);

        imguiDraw();
    }

    void onSound(AudioIOData &io) override { // Audio callback
        synthManager.render(io);
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
