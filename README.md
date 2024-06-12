This repository contains two things:
- the `spectrogram.hpp` file which can be included in any Allolib project and gives access to the class `al::Spectrogram`
- the `sound_design_with_spectrogram.cpp` file which demonstrates utilizing a Spectrogram object

# Spectrogram
The simplest way to use the Spectrogram object with an Allolib app is to put into any class that inherits `al::App`. In the class, it cannot be initialized outright and instead needs to be added as a member variable with an initializer list. 
The parameters needed in this list are the following:
- `int buffer_size`: determines how many past STFT outputs the Spetrogram should store and display
- `int fft_window_size`: determines how many bins the Short-Time Fourier Transform uses
- `int hop_size`: determines how many samples the STFT should wait for before calculating the next set of spectrum values
- `float x,y`: the x and y coordinate offset for the center of the spectrogram
- `float w,h`: the width and height of the spectrogram respectively
- `float ler`: legend ratio determines what percentage of the width should be used to display the spectrogram rather than the legend (values closer to 1 )
- `bool is_log`: boolean switch for linear or logarithmic frequency display
The values used in in the sound design file are a good example of useful settings for this.

For it to work properly, some member functions need to be called throughout the App class. These instructions assume the object has name `mySpec`:
- `onDraw(Graphics &g)`: call `mySpec.draw(g)` in this function to draw the Spectrogram object. When Allolib renders graphics, subsequent calls will overwrite previous ones. 
- `onSound(AudioIOData &io)`: there are two cases:
  - If there is a `while(io())` loop in the `onSound` function, call `mySpec.write_sample(io.out(0))` in it. The integer passed to the io object determines which audio channel the spectrogram will perform calculations on.
  - If there is no `while(io())` loop, and nothing else in `onSound` internally uses a while loop, call `mySpec.process_audio(io)` as the last function within `onSound`.

If these steps are followed correctly, the object will process and render the spectrogram for any input.

# Sound Design
