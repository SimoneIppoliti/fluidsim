#include "fluid.hpp"
#include "tools.hpp"
#include <cmath>
#include <algorithm>
#include <piksel/shape.hpp>
#include <glm/vec2.hpp>

#include <iostream>

namespace simulation {
    fluid::fluid() :
    size_(N), iter_(4), scale_(1), diff_(0.0f), visc_(0.0f), s_(N * N, 0.0f), density_(N * N, 0.0f),
    vX_(N * N, 0.0f), vX0_(N * N, 0.0f), vY_(N * N, 0.0f), vY0_(N * N, 0.0f)
    {}

    fluid::fluid(const int size, const int iter, const int scale, const float diff, const float visc) :
    size_(size), iter_(iter), scale_(scale), diff_(0.0f), visc_(0.0f), s_(size * size, 0.0f), density_(size * size, 0.0f),
    vX_(size * size, 0.0f), vX0_(size * size, 0.0f), vY_(size * size, 0.0f), vY0_(size * size, 0.0f)
    {}

    fluid::~fluid() {}

    void fluid::add_density(const int x, const int y, const float a) {
        density_[utils::IX(x, y, size_)] += a;
    }

    void fluid::reset_density(const int x, const int y) {
        density_[utils::IX(x, y, size_)] = 0.0f;
    }

    void fluid::add_velocity(const int x, const int y, const float aX, const float aY) {
        int i = utils::IX(x, y, size_);
        vX_[i] += aX;
        vY_[i] += aY;
    }

    void fluid::linear_solve(const int b, farr &x, farr &x0, const float a, const float c) {
        float c_recip = 1.0f / c;

        for (size_t t = 0; t < iter_; t++) {
            for (size_t j = 0; j < size_ - 1; j++) {
                for (size_t i = 0; i < size_ - 1; i++) {
                    x[utils::IX(i, j, size_)] = (
                        x0[utils::IX(i, j, size_)] +
                        a * (
                        x[utils::IX(i + 1, j, size_)] +
                        x[utils::IX(i - 1, j, size_)] +
                        x[utils::IX(i, j + 1, size_)] +
                        x[utils::IX(i, j - 1, size_)]
                        )
                    ) * c_recip;
                }
            }
            
            set_bounds(b, x);
        }
    }

    void fluid::project(farr &velX, farr &velY, farr &p, farr &div) {
        for (size_t j = 0; j < size_ - 1; j++) {
            for (size_t i = 0; i < size_ - 1; i++) {
                div[utils::IX(i, j, size_)] = -0.5f * (
                    velX[utils::IX(i + 1, j, size_)] -
                    velX[utils::IX(i - 1, j, size_)] +
                    velY[utils::IX(i, j + 1, size_)] -
                    velY[utils::IX(i, j - 1, size_)] ) / size_;
                p[utils::IX(i, j, size_)] = 0;
            }
        }
        
        set_bounds(0, div);
        set_bounds(0, p);
        linear_solve(0, p, div, 1.0f, 4.0f);

        for (size_t j = 0; j < size_ - 1; j++) {
            for (size_t i = 0; i < size_ - 1; i++) {
                velX[utils::IX(i, j, size_)] -= 0.5f * (p[utils::IX(i + 1, j, size_)] -
                p[utils::IX(i - 1, j, size_)]) * size_;

                velY[utils::IX(i, j, size_)] -= 0.5f * (p[utils::IX(i, j + 1, size_)] -
                p[utils::IX(i, j - 1, size_)]) * size_;
            }
        }

        set_bounds(1, velX);
        set_bounds(2, velY);
    }

    void fluid::advect(const int b, farr &d, farr &d0,  farr &velX, farr &velY, const float dt) {
        float i0, i1, j0, j1;
        
        float dtx = dt * (size_ - 2);
        float dty = dt * (size_ - 2);
        
        float s0, s1, t0, t1;
        float tmp1, tmp2, x, y;
        
        float Nfloat = size_;
        float jfloat, ifloat;
        size_t i, j;
        
        for(j = 1, jfloat = 1.0f; j < size_ - 1; j++, jfloat++) {
            for(i = 1, ifloat = 1.0f; i < size_ - 1; i++, ifloat++) {
                tmp1 = dtx * velX[utils::IX(i, j, size_)];
                tmp2 = dty * velY[utils::IX(i, j, size_)];
                x    = ifloat - tmp1;
                y    = jfloat - tmp2;
                
                if(x < 0.5f) x = 0.5f;
                if(x > Nfloat + 0.5f) x = Nfloat + 0.5f;
                i0 = std::floor(x);
                i1 = i0 + 1.0f;
                if(y < 0.5f) y = 0.5f;
                if(y > Nfloat + 0.5f) y = Nfloat + 0.5f;
                j0 = std::floor(y);
                j1 = j0 + 1.0f;
                
                s1 = x - i0;
                s0 = 1.0f - s1;
                t1 = y - j0;
                t0 = 1.0f - t1;
                
                int i0i = i0;
                int i1i = i1;
                int j0i = j0;
                int j1i = j1;
                
                d[utils::IX(i, j, size_)] =
                    s0 * (t0 * d0[utils::IX(i0i, j0i, size_)] + t1 * d0[utils::IX(i0i, j1i, size_)]) +
                    s1 * (t0 * d0[utils::IX(i1i, j0i, size_)] + t1 * d0[utils::IX(i1i, j1i, size_)]);
            }
        }
        set_bounds(b, d);
    }

    void fluid::diffuse(const int b, farr &x, farr &x0, const float diff, const float dt) {
        float a = dt * diff * (size_ - 2) * (size_ - 2);
        linear_solve(b, x, x0, a, 1 + 4 * a);
    }

    void fluid::fade_density(const float a) {
        for (size_t i = 0; i < density_.size(); i++) {
            float d = density_[i];
            density_[i] = std::max(d - a, 0.0f);
        }
    }

    void fluid::set_bounds(const int b, farr &x) {
        for(size_t i = 0; i < size_ - 1; i++) {
            x[utils::IX(i, 0, size_)] = b == 2 ? -x[utils::IX(i, 1, size_)] : x[utils::IX(i, 1, size_)];
            x[utils::IX(i, size_ - 1, size_)] = b == 2 ? -x[utils::IX(i, size_ - 2, size_)] : x[utils::IX(i, size_ - 2, size_)];
        }

        for(size_t j = 0; j < size_ - 1; j++) {
            x[utils::IX(0, j, size_)] = b == 1 ? -x[utils::IX(1, j, size_)] : x[utils::IX(1, j, size_)];
            x[utils::IX(size_ - 1, j, size_)] = b == 1 ? -x[utils::IX(size_ - 2, j, size_)] : x[utils::IX(size_ - 2, j, size_)];
        }
        
        x[utils::IX(0, 0, size_)] = 0.5f * (x[utils::IX(1, 0, size_)] + x[utils::IX(0, 1, size_)]);
        x[utils::IX(0, size_ - 1, size_)] = 0.5f * (x[utils::IX(1, size_ - 1, size_)] + x[utils::IX(0, size_ - 2, size_)]);
        x[utils::IX(size_ - 1, 0, size_)] = 0.5f * (x[utils::IX(size_ - 2, 0, size_)] + x[utils::IX(size_ - 1, 1, size_)]);
        x[utils::IX(size_ - 1, size_ - 1, size_)] = 0.5f * (x[utils::IX(size_ - 2, size_ - 1, size_)] + x[utils::IX(size_ - 1, size_ - 2, size_)]);
    }

    void fluid::step(const float dt) {
        diffuse(1, vX0_, vX_, visc_, dt);
        diffuse(2, vY0_, vY_, visc_, dt);
        
        project(vX0_, vY0_, vX_, vY_);
        
        advect(1, vX_, vX0_, vX0_, vY0_, dt);
        advect(2, vY_, vY0_, vX0_, vY0_, dt);
        
        project(vX_, vY_, vX0_, vY0_);
        
        diffuse(0, s_, density_, diff_, dt);
        advect(0, density_, s_, vX_, vY_, dt);
    }

    int fluid::size() const { return size_; }
    int fluid::scale() const { return scale_; }
    float fluid::density(const int x, const int y) const {
        return density_[utils::IX(x, y, size_)];
    }
    float fluid::velocity(const int x, const int y) const {
        float xx = vx(x, y);
        float yy = vy(x, y);
        return std::sqrt((xx * xx) + (yy * yy));
    }
    float fluid::vx(const int x, const int y) const {
        return std::abs(vX_[utils::IX(x, y, size_)] - vX0_[utils::IX(x, y, size_)]);
        //return std::abs(vX_[utils::IX(x, y, size_)]);
    }
    float fluid::vy(const int x, const int y) const {
        return std::abs(vY_[utils::IX(x, y, size_)] - vY0_[utils::IX(x, y, size_)]);
        //return std::abs(vY_[utils::IX(x, y, size_)]);
    }
}