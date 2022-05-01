#include "app.hpp"
#include "fluid.hpp"

int main() {
    int size = 64;
    int scale = 10;
    // size is 64, +2 accounts for broken borders
    simulation::fluid f(size + 2, 4, scale, 0.0f, 0.0f);

    App app(size * scale, size * scale, "Fluid simulation", false);

    app.set_fluid(f);
    app.start();
}