#ifndef IC_CAHNHILLARDMIXTURE_H_
#define IC_CAHNHILLARDMIXTURE_H_

#include "Set/Set.H"
#include "IC/IC.H"
#include "IO/ParmParse.H"

using namespace std;
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "mpi.h"
#include "IO/ParmParse.H"

namespace IC
{
class CahnHillardMixture : public IC
{
public:
    CahnHillardMixture(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}
    CahnHillardMixture(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : CahnHillardMixture(_geom)
    {
        pp.queryclass(name, *this);
    }

    void Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar)
    {
        Set::Vector DX(geom[lev].CellSize());
        Set::Vector size = Set::Size(geom[lev]);
        if (field_type == "phase") {
            for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx;
                bx = mfi.growntilebox();
                amrex::IndexType type = a_field[lev]->ixType();
                amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Vector x = Set::Position(i, j, k, geom[lev], type);
                    Set::Scalar min_grain_id = 0;
                    Set::Scalar d;
                    Set::Scalar d1 = std::numeric_limits<Set::Scalar>::infinity();

                    d = (x - X).lpNorm<2>();

                    if (geom[0].isPeriodic(0))
                    {
                        d1 = std::min((x - X + size(0) * Set::Vector::Unit(0)).lpNorm<2>(),
                            (x - X - size(0) * Set::Vector::Unit(0)).lpNorm<2>());
                    }
#if AMREX_SPACEDIM > 1
                    if (geom[0].isPeriodic(1))
                    {
                        d1 = std::min((x - X + size(1) * Set::Vector::Unit(1)).lpNorm<2>(),
                            (x - X - size(1) * Set::Vector::Unit(1)).lpNorm<2>());
                    }
#endif
#if AMREX_SPACEDIM > 2
                    if (geom[0].isPeriodic(2))
                    {
                        d1 = std::min((x - X + size(2) * Set::Vector::Unit(2)).lpNorm<2>(),
                            (x - X - size(2) * Set::Vector::Unit(2)).lpNorm<2>());
                    }
#endif
                    d = std::min(d, d1);
                    if (d <= (R + eps))
                    {
                        Set::Scalar m = 0.5 * (1 + erf((-d + R) / (eps)));
                        min_grain_id = min_grain_id + m * (1. - min_grain_id);
                    }

                    field(i, j, k, 1) = min_grain_id;
                    field(i, j, k, 2) = 0.98 - min_grain_id;

                    if (field(i, j, k, 2) < 0.) field(i, j, k, 2) = 0;
                });
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    field(i, j, k, 0) = 1. - field(i, j, k, 2) - field(i, j, k, 1);
                });
            }
        }
        else if (field_type == "concentration") {
            for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx;
                bx = mfi.growntilebox();
                amrex::IndexType type = a_field[lev]->ixType();
                amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Vector x = Set::Position(i, j, k, geom[lev], type);
                    Set::Scalar min_grain_id = 0;
                    Set::Scalar d;
                    Set::Scalar d1 = std::numeric_limits<Set::Scalar>::infinity();

                    d = (x - X).lpNorm<2>();

                    if (geom[0].isPeriodic(0))
                    {
                        d1 = std::min((x - X + size(0) * Set::Vector::Unit(0)).lpNorm<2>(),
                            (x - X - size(0) * Set::Vector::Unit(0)).lpNorm<2>());
                    }
#if AMREX_SPACEDIM > 1
                    if (geom[0].isPeriodic(1))
                    {
                        d1 = std::min((x - X + size(1) * Set::Vector::Unit(1)).lpNorm<2>(),
                            (x - X - size(1) * Set::Vector::Unit(1)).lpNorm<2>());
                    }
#endif
#if AMREX_SPACEDIM > 2
                    if (geom[0].isPeriodic(2))
                    {
                        d1 = std::min((x - X + size(2) * Set::Vector::Unit(2)).lpNorm<2>(),
                            (x - X - size(2) * Set::Vector::Unit(2)).lpNorm<2>());
                    }
#endif
                    d = std::min(d, d1);
                    if (d <= (R + eps))
                    {
                        Set::Scalar m = 0.5 * (1 + erf((-d + R) / (eps)));
                        min_grain_id = min_grain_id + m * (1. - min_grain_id);
                    }

                    field(i, j, k, 0) = min_grain_id;
                    field(i, j, k, 2) = 0.98 - min_grain_id;

                    if (field(i, j, k, 2) < 0.) field(i, j, k, 2) = 0;
                });
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    field(i, j, k, 1) = 1. - field(i, j, k, 0) - field(i, j, k, 2);
                });
            }
        }
        else Util::Abort(INFO, "Invalid field type");
        a_field[lev]->FillBoundary();
    }

private:
    std::string field_type; //seed
    std::string add_ic;

    Set::Vector X;
    Set::Scalar eps;
    Set::Scalar R; 
    
public:
    static void Parse(CahnHillardMixture& value, IO::ParmParse& pp)
    {
        pp.query_required("field.type", value.field_type); // Field type [concentration, parameter]
        pp_query_required("radius", value.R); // Radius of the sphere
        pp_queryarr("center", value.X); // Vector location of the sphere center
        pp_query_required("eps", value.eps); // Value of the field inside the sphere
    }
};
}
#endif
