#ifndef APP_HPP
#define APP_HPP

#include <piksel/baseapp.hpp>
#include <string>
//#include "fluid.hpp"

namespace simulation {
    class fluid;
}

class App : public piksel::BaseApp {
private:
    float dt_;
    float mx_, my_, pmx_, pmy_;

    bool ml_, mr_;

    simulation::fluid *fluid_;

public:
    App();
    App(int, int, std::string, bool);
    void setup();
    void draw(piksel::Graphics&);
    void mouseMoved(int, int);
    void mousePressed(int);
    void mouseReleased(int);

    void process_input();
    void set_fluid(simulation::fluid &);
    void frender(piksel::Graphics&) const;
};

#endif /* APP_HPP */