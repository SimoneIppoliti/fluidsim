#ifndef FLUID_HPP
#define FLUID_HPP

#include <vector>

namespace simulation {

    using farr = std::vector<float>;

    class fluid {
        protected:
            const int N = 256;
        
            int size_;
            int iter_;
            int scale_;

            float diff_;
            float visc_;

            std::vector<float> s_;
            std::vector<float> density_;

            std::vector<float> vX_;
            std::vector<float> vY_;

            std::vector<float> vX0_;
            std::vector<float> vY0_;

        public:
            fluid();
            fluid(const int, const int, const int, const float, const float);
            
            virtual ~fluid();

            void add_density(const int, const int, const float);
            void reset_density(const int, const int);
            void add_velocity(const int, const int, const float, const float);

            void linear_solve(const int, farr&, farr&, const float, const float);

            void project(farr&, farr&, farr&, farr&);
            void advect(const int, farr&, farr&,  farr&, farr&, const float);
            void diffuse(const int, farr&, farr&, const float, const float);
            void fade_density(const float);
            void set_bounds(const int, farr&);

            void step(const float);

            int size() const;
            int scale() const;
            float density(const int, const int) const;
            float velocity(const int, const int) const;
            float vx(const int, const int) const;
            float vy(const int, const int) const;
    };
}

#endif /* FLUID_HPP */