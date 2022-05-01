#include "app.hpp"
#include "fluid.hpp"

App::App() : piksel::BaseApp(640, 480), fluid_(nullptr),
dt_(0.0f), mx_(0.0f), my_(0.0f), pmx_(0.0f), pmy_(0.0f), ml_(false), mr_(false) {}

App::App(int width, int height, std::string title, bool fullscreen = false) :
piksel::BaseApp(width, height, title, fullscreen), fluid_(nullptr),
dt_(0.0f), mx_(0.0f), my_(0.0f), pmx_(0.0f), pmy_(0.0f), ml_(false), mr_(false) {}

void App::setup() {
    float now = millis();

    dt_ = millis() - now;
}

void App::draw(piksel::Graphics &g) {
    float now = millis();

    process_input();

    fluid_->step(dt_);

    g.background(glm::vec4(0.0f, 0.0f, 0.0f, 1.0f));
    frender(g);

    pmx_ = mx_;
    pmy_ = my_;
    dt_ = (millis() - now) / 1000.0f;
}

void App::mouseMoved(int x, int y) {
    mx_ = x; my_ = y;
}

void App::mousePressed(int button) {
    if (button == 0) ml_ = true;
    if (button == 1) mr_ = true;
}

void App::mouseReleased(int button) {
    if (button == 0) ml_ = false;
    if (button == 1) mr_ = false;
}

void App::process_input() {
    if (ml_) {
        int s = fluid_->scale();
        float x = mx_ / s + 1;
        float y = my_ / s + 1;
        fluid_->add_density(x, y, 250.0f * dt_);

        float vx = mx_ - pmx_;
        float vy = my_ - pmy_;
        fluid_->add_velocity(x, y, vx, vy);
    }
    if (mr_) {
        int s = fluid_->scale();
        float x = mx_ / s + 1;
        float y = my_ / s + 1;
        fluid_->reset_density(x, y);
    }
}

void App::set_fluid(simulation::fluid &f) { fluid_ = &f; }

void App::frender(piksel::Graphics& g) const {
    int s = fluid_->size();
    int t = fluid_->scale();

    g.noStroke();

    for (size_t j = 0; j < s; j++) {
        for (size_t i = 0; i < s; i++) {
            int x = i * t;
            int y = j * t;
            float d = fluid_->density(i, j);
            g.fill(glm::vec4(d, 0.1f, 100.0f / d, d));
            g.rect(x - t, y - t, t, t); // -t accounts for broken borders
        }
    }

    fluid_->fade_density(0.3f * dt_);
}