#include <eigen3/Eigen/Core>
#include "Set/Set.H"

/// A bunch of solvers
namespace Solver
{
/// Local solvers
namespace Local
{

namespace Riemann
{

/// HLLC Riemann Solver
class HLLC
{
public:
    struct State {
        Set::Scalar rho;
        Set::Scalar velocity_normal;
        Set::Scalar velocity_tangent;
        Set::Scalar pressure;
        // Construtor for convenience
        State()
        {
            rho = 0.0;
            velocity_normal = 0.0;
            velocity_tangent = 0.0;
            pressure = 0.0;
        }
        State(Set::Scalar a_rho, Set::Scalar a_velocity_normal, Set::Scalar a_velocity_tangent, Set::Scalar a_pressure)
            : rho(a_rho),
            velocity_normal(a_velocity_normal), velocity_tangent(a_velocity_tangent),
            pressure(a_pressure) {}
    };

    struct Flux {
        Set::Scalar Mass;
        Set::Scalar Energy;
        Set::Scalar Momentum_normal;
        Set::Scalar Momentum_tangent;
    };


    // todo: convert to static functor by replacing with static () when c++23 is released
    static Flux Solve(State lo, State hi, Set::Scalar gamma) {

        //lo primitive variables
        Set::Scalar uNorm_lo = lo.velocity_normal;
        Set::Scalar uTang_lo = lo.velocity_tangent;
        Set::Scalar ke_lo = 0.5 * lo.rho * (uNorm_lo * uNorm_lo + uTang_lo * uTang_lo);
        Set::Scalar ue_lo = lo.pressure / (gamma - 1);

        //hi primitive variables
        Set::Scalar uNorm_hi = hi.velocity_normal;
        Set::Scalar uTang_hi = hi.velocity_tangent;
        Set::Scalar ke_hi = 0.5 * hi.rho * (uNorm_hi * uNorm_hi + uTang_hi * uTang_hi);
        Set::Scalar ue_hi = hi.pressure / (gamma - 1);

        //Find the largest eigenvalues in the normal direction to the interface
        Set::Scalar a_lo = std::sqrt(gamma * lo.pressure / lo.rho);
        Set::Scalar a_hi = std::sqrt(gamma * hi.pressure / hi.rho);

        //Compute HLL wave speed
        Set::Scalar S_min = std::min(uNorm_lo, uNorm_hi) - max(a_lo, a_hi);
        Set::Scalar S_max = std::max(uNorm_lo, uNorm_hi) + max(a_lo, a_hi);

        //Compute lagrangian sound speed
        Set::Scalar aLag_min = lo.rho * (uNorm_lo - S_min);
        Set::Scalar aLag_max = hi.rho * (S_max - uNorm_hi);

        //Compute acoustic star state
        Set::Scalar ustar = (aLag_max * uNorm_hi + aLag_min * uNorm_lo + (lo.pressure - hi.pressure)) / (aLag_max + aLag_min);
        Set::Scalar pstar = (aLag_max * lo.pressure + aLag_min * hi.pressure + aLag_min * aLag_max * (uNorm_lo - uNorm_hi)) / (aLag_max + aLag_min);

        //Left star region variables
        Set::Scalar rstar_lo = lo.rho * (S_min - uNorm_lo) / (S_min - ustar);
        Set::Scalar uestar_lo = ((S_min - uNorm_lo) * ue_lo - lo.pressure * uNorm_lo + pstar * ustar) / (S_min - ustar);
        Set::Scalar kestar_lo = ((S_min - uNorm_lo) * ke_lo - lo.pressure * uNorm_lo + pstar * ustar) / (S_min - ustar);

        //Right star region variables
        Set::Scalar rstar_hi = hi.rho * (S_max - uNorm_hi) / (S_max - ustar);
        Set::Scalar uestar_hi = ((S_max - uNorm_hi) * ue_hi - hi.pressure * uNorm_hi + pstar * ustar) / (S_max - ustar);
        Set::Scalar kestar_hi = ((S_max - uNorm_hi) * ke_hi - hi.pressure * uNorm_hi + pstar * ustar) / (S_max - ustar);

        //Sample the solution at x/t=0
        Set::Scalar rho_face, u_face, p_face, ke_face, ue_face;
        if (S_min > 0.0) {
            rho_face = lo.rho;
            u_face = uNorm_lo;
            p_face = lo.pressure;
            ue_face = ue_lo;
            ke_face = ke_lo;
        }
        else if (ustar > 0.0) {
            rho_face = rstar_lo;
            u_face = ustar;
            p_face = pstar;
            ue_face = uestar_lo;
            ke_face = kestar_lo;
        }
        else if (S_max > 0.0) {
            rho_face = rstar_hi;
            u_face = ustar;
            p_face = pstar;
            ue_face = uestar_hi;
            ke_face = kestar_hi;
        }
        else {
            rho_face = hi.rho;
            u_face = uNorm_hi;
            p_face = hi.pressure;
            ue_face = ue_hi;
            ke_face = ke_hi;
        };

        //Compute the Godunov flux
        Flux fl;
        fl.Mass = rho_face * u_face;
        fl.Momentum_normal = (rho_face * u_face * u_face) + p_face;
        if (fl.Mass > 0.0) {
            fl.Momentum_tangent = (fl.Mass * uTang_lo) * u_face;
        }
        else {
            fl.Momentum_tangent = (fl.Mass * uTang_hi) * u_face;
        };
        fl.Energy = (ke_face + ue_face) * u_face;

        return fl;

    }



    static int Test()
    {
        int failed = 0;
        for (int i = 0; i < 1; i++)
        {

            State hi, lo, fixed, pos_vel1, pos_vel2, neg_vel1, neg_vel2, comp1, comp2;
            hi.rho = Util::Random();
            hi.velocity_normal = Util::Random();
            hi.velocity_tangent = Util::Random();
            hi.pressure = Util::Random();

            lo.rho = Util::Random();
            lo.velocity_normal = Util::Random();
            lo.velocity_tangent = Util::Random();
            lo.pressure = Util::Random();

            fixed.rho = Util::Random();
            fixed.velocity_normal = 0.0f;
            fixed.velocity_tangent = Util::Random();
            fixed.pressure = Util::Random();

            pos_vel1.rho = Util::Random();
            pos_vel1.velocity_normal = Util::Random();
            pos_vel1.velocity_tangent = Util::Random();
            pos_vel1.pressure = Util::Random();

            neg_vel1.rho = pos_vel1.rho;
            neg_vel1.velocity_normal = -pos_vel1.velocity_normal;
            neg_vel1.velocity_tangent = -pos_vel1.velocity_tangent;
            neg_vel1.pressure = pos_vel1.pressure;

            pos_vel2.rho = Util::Random();
            pos_vel2.velocity_normal = Util::Random();
            pos_vel2.velocity_tangent = Util::Random();
            pos_vel2.pressure = Util::Random();

            neg_vel2.rho = pos_vel2.rho;
            neg_vel2.velocity_normal = -pos_vel2.velocity_normal;
            neg_vel2.velocity_tangent = -pos_vel2.velocity_tangent;
            neg_vel2.pressure = pos_vel2.pressure;

            comp1.rho = 1.0f;
            comp1.velocity_normal = 2.0f;
            comp1.velocity_tangent = 3.0f;
            comp1.pressure = 4.0f;

            comp2.rho = 1.1f;
            comp2.velocity_normal = 2.2f;
            comp2.velocity_tangent = 3.3f;
            comp2.pressure = 4.4f;


            // 
            // Zero flux test: if the state is the same on both sides,
            //                 the flux should be zero.
            // 
            Flux zero = Solve(fixed, fixed, 1.4);
            if (fabs(zero.Mass) > 1E-10) { Util::Warning(INFO, "Nonzero mass flux ", zero.Mass); failed++; }
            if (fabs(zero.Momentum_normal) > 1E-10) { Util::Warning(INFO, "Nonzero x momentum flux ", zero.Momentum_normal); failed++; }
            if (fabs(zero.Momentum_tangent) > 1E-10) { Util::Warning(INFO, "Nonzero y momentum flux ", zero.Momentum_tangent); failed++; }
            if (fabs(zero.Energy) > 1E-10) { Util::Warning(INFO, "Nonzero energy flux ", zero.Energy); failed++; }

            // 
            // Antisymmetry test: reversing hi and lo should change the sign
            // of the flux.
            // 

            Flux plus = Solve(hi, lo, 1.4);
            Flux minus = Solve(lo, hi, 1.4);

            if (fabs(plus.Mass + minus.Mass) > 1E-10)
            {                
                Util::Warning(INFO, "mass flux error: plus=", plus.Mass, " minus=", minus.Mass); failed++;            
            }
            if (fabs(plus.Momentum_normal + minus.Momentum_normal) > 1E-10)
            {
                Util::Warning(INFO, "momentum x flux error: plus=", plus.Momentum_normal, " minus=", minus.Momentum_normal); failed++;
            }
            if (fabs(plus.Momentum_tangent + minus.Momentum_tangent) > 1E-10)
            {
                Util::Warning(INFO, "momentum y flux error: plus=", plus.Momentum_tangent, " minus=", minus.Momentum_tangent); failed++;
            }
            if (fabs(plus.Energy + minus.Energy) > 1E-10)
            {
                Util::Warning(INFO, "energy flux error: plus=", plus.Energy, " minus=", minus.Energy); failed++;
            }


            // 
            // Wave Propogation Test: changing the sign of velocity should produce negative mass
            //                   and energy fluxes and the same momentum flux
            // 

            Flux pos = Solve(pos_vel1, pos_vel2, 1.4);
            Flux neg = Solve(neg_vel1, neg_vel2, 1.4);

            if (fabs(pos.Mass + neg.Mass) > 1E-10)
            {
                Util::Warning(INFO, "mass flux error: positive vel=", pos.Mass, " negative vel=", neg.Mass); failed++;
            }
            if (fabs(pos.Momentum_normal - neg.Momentum_normal) > 1E-10)
            {
                Util::Warning(INFO, "momentum x flux error: positive vel=", pos.Momentum_normal, " negative vel=", neg.Momentum_normal); failed++;            
            }
            if (fabs(pos.Momentum_tangent - neg.Momentum_tangent) > 1E-10)
            {
                Util::Warning(INFO, "momentum y flux error: positive vel=", pos.Momentum_tangent, " negative vel=", neg.Momentum_tangent); failed++;
            }
            if (fabs(pos.Energy + neg.Energy) > 1E-10)
            {
                Util::Warning(INFO, "energy flux error: positive vel=", pos.Energy, " negative vel=", neg.Energy); failed++;
            }

            Flux comp = Solve(comp1, comp2, 1.4);
        }

        return failed;
    }

};
}
}
}
