#include "ChemReduction.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "Numeric/Function.H"
#include <cmath>

namespace Integrator
{

ChemReduction::ChemReduction() {}

ChemReduction::ChemReduction(IO::ParmParse& pp) : ChemReduction()
{
    pp.queryclass(*this);
}

void
ChemReduction::Parse(ChemReduction& value, IO::ParmParse& pp)
{
    BL_PROFILE("Integrator::ChemReduction::ChemReduction()");
    {
        { // Parse Properties
            pp.query_required("molar_volume", value.V); // Domain molar volume 
            pp.query_required("temperature", value.T); // Domain temperature (Isothermal)
            pp.query_default("gas_constant", value.R, 8.3145); // Universal gas constant
            pp.query_required("reaction_constant", value.RC); // Reaction Constant

            pp.query_default("refinement_criterion", value.phi_refinement_criterion, 0.001); // Criterion for mesh refinement
            pp.query_default("ghost_cells_count", value.ghost_count, 2); // Number of ghost cells in each field

            pp.query_required("wustite.bulk", value.wustite.K); // Wustite bulk modulus 
            pp.query_required("wustite.shear", value.wustite.G); // Wustite shear modulus
            pp.query_required("wustite.mobility", value.wustite.L); // Wustite mobility coeficient
            pp.query_required("wustite.diffusivity.oxygen", value.wustite.D.O); // Wustite oxygen diffusion coeficient
            pp.query_required("wustite.diffusivity.hydrogen", value.wustite.D.H); // Wustite hydrogen diffusion coeficient
            pp.query_required("wustite.diffusivity.water", value.wustite.D.W); // Wustite water diffusion coeficient 
            pp.query_required("wustite.energy.gradient", value.wustite.E.G); // Wustite energy gradient
            pp.query_required("wustite.energy.barrier", value.wustite.E.B); // Wustite energy barrier

            pp.query_required("ferrite.bulk", value.ferrite.K); // Ferrite bulk modulus 
            pp.query_required("ferrite.shear", value.ferrite.G); // Ferrite shear modulus 
            pp.query_required("ferrite.mobility", value.ferrite.L); // Ferrite mobility coeficient
            pp.query_required("ferrite.diffusivity.oxygen", value.ferrite.D.O); // Ferrite oxygen diffusion coeficient
            pp.query_required("ferrite.diffusivity.hydrogen", value.ferrite.D.H);// Ferrite hydrogen diffusion coeficient
            pp.query_required("ferrite.diffusivity.water", value.ferrite.D.W);// Ferrite water diffusion coeficient
            pp.query_required("ferrite.energy.gradient", value.ferrite.E.G);// Ferrite energy gradient
            pp.query_required("ferrite.energy.barrier", value.ferrite.E.B);// Ferrite energy barrier

            pp.query_required("gas.bulk", value.gas.K); // Water/gas bulk modulus
            pp.query_required("gas.shear", value.gas.G); // Water/gas shear modulus
            pp.query_required("gas.mobility", value.gas.L); // Water/gas mobility coeficient 
            pp.query_required("gas.diffusivity.oxygen", value.gas.D.O); // Water/gas oxygen diffusion coeficient
            pp.query_required("gas.diffusivity.hydrogen", value.gas.D.H); // Water/gas hydrogen diffusion coeficient
            pp.query_required("gas.diffusivity.water", value.gas.D.W); // Water/gas water diffusion coeficient 
            pp.query_required("gas.energy.gradient", value.gas.E.G); // Water/gas energy gradient
            pp.query_required("gas.energy.barrier", value.gas.E.B); // Water/gas energy barrier
        } // End Parse Properties

        { // Allen-Cahn (Non-Conserved) fields
            value.bc_phi = new BC::Constant(3);
            pp.queryclass("pf.phi.bc", *static_cast<BC::Constant*>(value.bc_phi)); // See :ref:`BC::Constant`
            value.RegisterNewFab(value.phi_mf, value.bc_phi, 3, value.ghost_count, "phi", true);
            value.RegisterNewFab(value.phi_o_mf, value.bc_phi, 3, value.ghost_count, "phi_old", false);
            value.ic_phi = new IC::CahnHillardMixture(value.geom, pp, "phi.ic.cahnhillard");
        } // End Allen-Cahn Register Fab

        { // Cahn-Hillard (Conserved) fields
            value.bc_chi = new BC::Constant(3);
            pp.queryclass("pf.chi.bc", *static_cast<BC::Constant*>(value.bc_chi)); // See :ref:`BC::Constant`
            value.RegisterNewFab(value.chi_mf, value.bc_chi, 3, value.ghost_count, "chi", true);
            value.RegisterNewFab(value.chi_o_mf, value.bc_chi, 3, value.ghost_count, "chi_old", false);
            value.ic_chi = new IC::CahnHillardMixture(value.geom, pp, "chi.ic.cahnhillard");
        } // End Cahn-Hillard Register Fab       
    } // End BL_PROFILE
} // End Parse Function 

void ChemReduction::Initialize(int lev)
{
    BL_PROFILE("Integrator::ChemReduction::Initialize");
    {
        ferrite.M.O = (ferrite.D.O * V) / (R * T);
        ferrite.M.H = (ferrite.D.H * V) / (R * T);
        ferrite.M.W = (ferrite.D.W * V) / (R * T);

        wustite.M.O = (wustite.D.O * V) / (R * T);
        wustite.M.H = (wustite.D.H * V) / (R * T);
        wustite.M.W = (wustite.D.W * V) / (R * T);

        gas.M.O = (gas.D.O * V) / (R * T);
        gas.M.H = (gas.D.H * V) / (R * T);
        gas.M.W = (gas.D.W * V) / (R * T);

        ic_phi->Initialize(lev, phi_mf);
        ic_chi->Initialize(lev, chi_o_mf);

        ic_chi->Initialize(lev, chi_mf);
        ic_chi->Initialize(lev, chi_o_mf);

    } // End BL_PROFILE
} // End Initialize Function

void ChemReduction::Advance(int lev, Set::Scalar, Set::Scalar dt)
{
    BL_PROFILE("Integrador::ChemReduction::Advance");
    {
        const Set::Scalar* DX = geom[lev].CellSize();

        std::swap(phi_mf[lev], phi_o_mf[lev]);
        std::swap(chi_mf[lev], chi_o_mf[lev]);

        // Free Energy Functions
        Numeric::Function::Polynomial<4> f_f(-40372, -36.031, 90.573, -182.51, 163.708); // wustite free energy
        Numeric::Function::Polynomial<2> f_w(41.953, -213.46, 173.263); // ferrite free energy
        // Interpolation Functions
        Numeric::Function::Polynomial<5> p(0., 0., 0., 10., -15., 6.); // interpolation function
        Numeric::Function::Polynomial<4> g(0., 0., 1., -2., 1.); // interpolation function

        Numeric::Function::Polynomial<4> dp(0., 0., 10.*2., -15.*3., 4.*6.);
        Numeric::Function::Polynomial<3> dg(0., 1., -2.*2., 1.*3.);

        for (amrex::MFIter mfi(*phi_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();

            amrex::Array4<Set::Scalar> const& phi = (*phi_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& phio = (*phi_o_mf[lev]).array(mfi);

            amrex::Array4<Set::Scalar> const& chi = (*chi_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& chio = (*chi_o_mf[lev]).array(mfi);

            // Parallel Loop for Cahn-Hillard Fields
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar lap_chiH = Numeric::Laplacian(chio, ferrite.D.H, i, j, k, 0, DX);
                Set::Scalar lap_chiW = Numeric::Laplacian(chio, ferrite.D.W, i, j, k, 1, DX);
                Set::Scalar lap_chiO = Numeric::Laplacian(chio, ferrite.D.O, i, j, k, 2, DX);

                Set::Scalar kfor = g(phi(i, j, k, 2)) * RC * chio(i, j, k, 2) * chio(i, j, k, 0);

                chi(i, j, k, 0) = chio(i, j, k, 0) - dt * (lap_chiH - 2 * kfor); // Hydrogen
                chi(i, j, k, 1) = chio(i, j, k, 1) - dt * (lap_chiO - kfor); // Water
                chi(i, j, k, 2) = chio(i, j, k, 2) - dt * (lap_chiW + kfor); // Oxygen
            });

            // Parallel Loop for Allen-Cahn Fields
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar df_divf = Numeric::Laplacian(phi, ferrite.E.G, i, j, k, 0, DX) - Numeric::Laplacian(phi, gas.E.G, i, j, k, 2, DX);
                    Set::Scalar df_divw = Numeric::Laplacian(phi, wustite.E.G, i, j, k, 1, DX) - Numeric::Laplacian(phi, gas.E.G, i, j, k, 2, DX);

                    Set::Scalar eta = p(phio(i, j, k, 0)) * f_f(chio(i, j, k, 2)) + p(phi(i, j, k, 2)) * f_w(chio(i, j, k, 2)); // eta_h(chiH) + eta_w(chiW)
                    Set::Scalar df_dphif = eta * (dp(phio(i, j, k, 0)) - dp(phi(i, j, k, 2))) + dg(phio(i, j, k, 0)) * ferrite.E.B - dg(phi(i, j, k, 2)) * gas.E.B;
                    Set::Scalar df_dphiw = eta * (dp(phio(i, j, k, 1)) - dp(phi(i, j, k, 2))) + dg(phio(i, j, k, 1)) * ferrite.E.B - dg(phi(i, j, k, 2)) * gas.E.B;

                    phi(i, j, k, 0) = phio(i, j, k, 0) + dt * ferrite.L * (df_divf - df_dphif); // Ferrite
                    phi(i, j, k, 1) = phio(i, j, k, 1) + dt * wustite.L * (df_divw - df_dphiw); // Wustite
                });

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    phi(i, j, k, 2) = 1.0 - phi(i, j, k, 0) - phi(i, j, k, 1); // Gas
                });
            }
        }// End For Loop
    } // End BL_PROFILE
} //Function

void ChemReduction::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar, int )
{
    BL_PROFILE("Integrator::ChemReduction::TagCellsForRefinement");
    {
        const Set::Scalar* DX = geom[lev].CellSize();
        Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

        // Phi criterion for refinement
        for (amrex::MFIter mfi(*phi_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<char> const& tags = a_tags.array(mfi);
            amrex::Array4<const Set::Scalar> const& phi = (*phi_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector gradf = Numeric::Gradient(phi, i, j, k, 0, DX);
                Set::Vector gradw = Numeric::Gradient(phi, i, j, k, 1, DX);
                if (gradf.lpNorm<2>() * dr * 2 > phi_refinement_criterion || gradw.lpNorm<2>() * dr * 2 > phi_refinement_criterion)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    } // End BL_PROFILE
} // End TagCell

void ChemReduction::Regrid(int lev, Set::Scalar)
{
    BL_PROFILE("Integrator::ChemReduction::Regrid");
    {
        if (lev < finest_level) return;
        Util::Message(INFO, "Regridding on level", lev);
    } // End BL_PROFILE
} // End Regrid 

} // namespace Integrator


