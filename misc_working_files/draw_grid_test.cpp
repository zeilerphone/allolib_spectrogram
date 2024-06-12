#include "al/app/al_App.hpp"
#include "Gamma/Oscillator.h"

#include "al/graphics/al_Shapes.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"

#include <cstdio>

using namespace al;
using namespace std;


class Strip : public Mesh {
public:
    Strip(float w, float h, int s, RGB cStart, RGB cEnd) :
        width(w), height(h), segments(s), segment_height(h / s),
        color_start(cStart), color_end(cEnd), color_dt((cEnd - cStart) / s)
    {
        this->primitive(TRIANGLE_STRIP);
        this->vertex(-width / 2, -height / 2);
        this->vertex(width / 2, -height / 2);
        for (int i = 0; i <= segments; i++) {
            this->color(color_start);
            this->vertex(-width / 2, -height / 2 + segment_height * static_cast<float>(i));
            this->color(color_start);
            this->vertex(width / 2, -height / 2 + segment_height * static_cast<float>(i));
        }
        this->color(color_start);
        this->vertex(-width / 2, height / 2);
        this->color(color_start);
        this->vertex(width / 2, height / 2);
    }

    void map_value_to_color(vector<float> vals) {
        // need to pass in vector of values between zero and one
        // zero is for no amplitude, one is for max
        // first value is the lowest pitch bucket, last is highest
        // vector needs to have size 'segments' or 'segments + 1'
        for (int i = 0; i <= segments; i++) {
            RGB curCol = color_start + (color_end - color_start) * vals[i];
            this->colors()[2 * i] = curCol;
            this->colors()[2 * i + 1] = curCol;
        }
    }

protected:
    float width, height, segment_height;
    int segments;
    RGB color_start, color_end, color_dt;
};


class Grid {
public:
    Grid(float w, float h, int ws, int hs, RGB low, RGB high) :
        width(w), height(h), width_of_strips(w / ws),
        width_segments(ws), height_segments(hs),
        color_low(low), color_high(high), color_delta(high-low / hs)
    {
        for (int i = 0; i < width_segments; i++) {
            Strip thisStrip(width_of_strips, height, height_segments, color_low, color_high);
            strips.push_back(thisStrip);
        }
        writePointer = 0;
    }

    /*~Grid() {
        for (int i = 0; i < width_segments; i++) {
            delete strips[i];
        }
    }*/

    void write_data(vector<float> vals) {
        strips[writePointer].map_value_to_color(vals);
        writePointer = (writePointer + 1) % width_segments;
    }

    Strip read_data(int index) {
        return strips[(writePointer + index) % width_segments];
    }


protected:
    float width, height, width_of_strips;
    int width_segments, height_segments;
    RGB color_low, color_high, color_delta;
    int writePointer;
    vector<Strip> strips;
};

//class Pixel : public Mesh {
//public:
//    Pixel(float x, float y, float w, float h, RGB c) :
//        xPos(x), yPos(y), width(w), height(h), color(c)
//    {
//        this->primitive(Triangle_Strip);
//        this->vertex(xPos + (-width / 2), yPos + ( height / 2));
//        this->vertex(xPos + ( width / 2), yPos + ( height / 2));
//        this->vertex(xPos + ( width / 2), yPos + (-height / 2));
//        this->vertex(xPos + (-width / 2), yPos + (-height / 2));
//        this->color(RGB);
//    }
//
//protected:
//    float xPos, yPos, width, height;
//    RGB color;
//};

//class Grid2 : public Mesh {
//public:
//    Grid2(float w, float h, int ws, int hs) :
//        width(w), height(h), width_D(ws), height_D(hs)
//    {
//         
//    }
//
//protected:
//    float width, height;
//    int width_D, height_D;
//};

struct MyApp: public App {
    RGB starting_color{ 0.1f, 0.1f, 0.1f };
    RGB ending_color{ 0.2f, 0.8f, 0.5f };
    gam::Sine<> oscillator;
    int number_of_strips_across = 200;
    int number_of_strips_high = 128;
    float width_of = 1.9;
    Grid grid{ width_of, 1.9, number_of_strips_across, number_of_strips_high, starting_color, ending_color };
    // Strip strip{ 0.01, 1.9, 2048, starting_color, ending_color};
    vector<float> tempVals;

    void onInit() override { // Called on app start
        tempVals.resize(number_of_strips_high + 1);
        std::cout << "onInit()" << std::endl;
        for (int i = 0; i <= number_of_strips_high; i++) {
            tempVals[i] = (static_cast<float>(i) / static_cast<float>(number_of_strips_across));
        }
    }

    void onCreate() override { // Called when graphics context is available
        std::cout << "onCreate()" << std::endl;
    }

    void onAnimate(double dt) override { // Called once before drawing
        //color += 0.01;
        //if (color > 1.0) {
        //    color -= 1.0;
        //}
        oscillator.freq(220);
        for (int i = 0; i <= number_of_strips_high; i++) {
            tempVals[i] += 1.f / number_of_strips_high;
            if (tempVals[i] >= 1) tempVals[i] = 0;
        }
        grid.write_data(tempVals);
    } 

    void onDraw(Graphics &g) override { // Draw function
        //g.clear(color);
        g.clear();
        g.color(1);
        g.camera(Viewpoint::IDENTITY);
        g.meshColor();
        g.translate((-width_of + width_of / number_of_strips_across) / 2.f, 0);
        for (int i = 0; i < number_of_strips_across; i++) {
            g.draw(grid.read_data(i));
            g.translate(width_of / number_of_strips_across, 0);
        }
    }

    void onSound(AudioIOData &io) override { // Audio callback
        while (io()) {
            io.out(0) = oscillator() * 0.1;
        }
    }

    void onMessage(osc::Message &m) override { // OSC message callback
        m.print();
    }
};


int main()
{
    MyApp app;
    auto dev = AudioDevice("devicename");
    app.configureAudio(dev, 44100, 256, 2, 2);

    app.start();
    return 0;
}
