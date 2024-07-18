#include "Hydro.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "BC/Nothing.H"
#include "Numeric/Stencil.H"
#include "Numeric/Function.H"
#include "IC/Laminate.H"
#include "IC/PSRead.H"
#include "IC/Expression.H"
#include "IC/BMP.H"
#include "IC/PNG.H"
#include "Solver/Local/Riemann/Roe.H"

namespace Integrator
{

Hydro::Hydro(IO::ParmParse& pp) : Hydro()
{
    pp_queryclass(*this);
}

void
Hydro::Parse(Hydro& value, IO::ParmParse& pp)
{
    BL_PROFILE("Integrator::Hydro::Hydro()");
    {
        // unused...?
        pp_query_required("r_refinement_criterion",     value.r_refinement_criterion    );
        // refinement based on energy flux
        pp_query_required("e_refinement_criterion",     value.e_refinement_criterion    );
        // refinement based on mass flux
        pp_query_required("m_refinement_criterion",     value.m_refinement_criterion    );
        // refinement based on :math:`\eta`
        pp_query_required("eta_refinement_criterion",   value.eta_refinement_criterion  );
        // refinement based on vorticity (disabled by default)
        pp_query_default("omega_refinement_criterion", value.omega_refinement_criterion, 1E100);

        // ideal gas :math:`\gamma`
        pp_query_required("gamma", value.gamma);
        // courant-friedrichs-levy ncondition
        pp_query_required("cfl", value.cfl);
        // viscosity
        pp_query_required("mu", value.mu);

        pp_forbid("rho.bc","replaced with density.bc");
        value.density_bc = new BC::Constant(1, pp, "density.bc");
        pp_forbid("p.bc","replaced with pressure.bc");
        value.pressure_bc = new BC::Constant(1, pp, "pressure.bc");
        pp_forbid("v.bc","replaced with velocity.bc");
        value.velocity_bc = new BC::Constant(2, pp, "velocity.bc");
        value.eta_bc = new BC::Constant(1, pp, "pf.eta.bc");

        value.momentum_bc = new BC::Constant(2, pp, "momentum.bc");
        value.energy_bc = new BC::Constant(1, pp, "energy.bc");


    }
    // Register FabFields:
    {
        int nghost = 2;

        value.RegisterNewFab(value.eta_mf,           value.eta_bc,      1, nghost, "eta",          true );
        value.RegisterNewFab(value.eta_old_mf,       value.eta_bc,      1, nghost, "eta_old",      false);
        value.RegisterNewFab(value.etadot_mf,        value.eta_bc,      1, nghost, "etadot",       true );
        value.RegisterNewFab(value.density_mf,       value.density_bc,  1, nghost, "density",      true );
        value.RegisterNewFab(value.density_old_mf,   value.density_bc,  1, nghost, "rho_old",      false);
        value.RegisterNewFab(value.energy_mf,        value.energy_bc,   1, nghost, "energy",       true );
        value.RegisterNewFab(value.energy_old_mf,    value.energy_bc,   1, nghost, "energy_old" ,  false);
        value.RegisterNewFab(value.momentum_mf,      value.momentum_bc, 2, nghost, "momentum",     true );
        value.RegisterNewFab(value.momentum_old_mf,  value.momentum_bc, 2, nghost, "momentum_old", false);
        value.RegisterNewFab(value.velocity_mf,      value.velocity_bc, 2, nghost, "velocity",     true );
        value.RegisterNewFab(value.velocity_old_mf,  value.velocity_bc, 2, nghost, "velocity_old", false);
        value.RegisterNewFab(value.vorticity_mf,    &value.bc_nothing,  1, nghost, "vorticity",    true );
        value.RegisterNewFab(value.pressure_mf,      value.pressure_bc, 1, nghost, "pressure",     true );
        value.RegisterNewFab(value.pressure_old_mf,  value.pressure_bc, 1, nghost, "pressure_old",     false);
        value.RegisterNewFab(value.vInjected_mf,    &value.bc_nothing,  2, nghost, "vInjected",    true );
        value.RegisterNewFab(value.rhoInterface_mf, &value.bc_nothing,  1, nghost, "rhoInterface", true );
        value.RegisterNewFab(value.q_mf,            &value.bc_nothing,  2, nghost, "q",            true );

        value.RegisterNewFab(value.solid.velocity_mf, &value.bc_nothing, 2, nghost, "solid.velocity",  true);
        value.RegisterNewFab(value.solid.density_mf,  &value.bc_nothing, 1, nghost, "solid.density",  true);
        value.RegisterNewFab(value.solid.pressure_mf,  &value.bc_nothing, 1, nghost, "solid.pressure",  true);
        value.RegisterNewFab(value.solid.momentum_mf, &value.bc_nothing, 2, nghost, "solid.momentum", true);
        value.RegisterNewFab(value.solid.energy_mf,   &value.bc_nothing, 1, nghost, "solid.energy",   true);

        value.RegisterNewFab(value.source_mf, &value.bc_nothing, 4, nghost, "source", true);
    }
    {
        std::string type = "constant";
        // initial condition for :math:`\eta`
        pp_query_validate("eta.ic.type", type,{"constant","laminate","expression","bmp","png"});
        if (type == "constant") value.eta_ic = new IC::Constant(value.geom, pp, "eta.ic.constant");
        else if (type == "laminate") value.eta_ic = new IC::Laminate(value.geom, pp, "eta.ic.laminate");
        else if (type == "expression") value.eta_ic = new IC::Expression(value.geom, pp, "eta.ic.expression");
        else if (type == "bmp") value.eta_ic = new IC::BMP(value.geom, pp, "eta.ic.bmp");
        else if (type == "png") value.eta_ic = new IC::PNG(value.geom, pp, "eta.ic.png");
    }
    {
        std::string type;
        pp_forbid("Velocity.ic","replaced with velocity.ic");
        // initial condition for velocity
        pp_query_validate("velocity.ic.type", type,{"constant","expression"});
        if (type == "constant") value.velocity_ic = new IC::Constant(value.geom, pp, "velocity.ic.constant");
        else if (type == "expression") value.velocity_ic = new IC::Expression(value.geom, pp, "velocity.ic.expression");
        else Util::Abort(INFO, "Invalid Velocity.ic: ", type);
    }
    {
        std::string type = "constant";
        pp_forbid("Pressure.ic","replaced with pressure.ic");
        // initial condition for pressure
        pp_query_validate("pressure.ic.type", type, {"constant","expression"});
        if (type == "constant") value.pressure_ic = new IC::Constant(value.geom, pp, "pressure.ic.constant");
        else if (type == "expression") value.pressure_ic = new IC::Expression(value.geom, pp, "pressure.ic.expression");
        else Util::Abort(INFO, "Invalid Pressure.ic: ", type);
    }
    {
        std::string type = "constant";
        pp_forbid("SolidVelocity.ic","replaced by solid.velocity");
        // initial condition for velocity in the solid phase
        pp_query_validate("solid.velocity.ic.type", type, {"constant","expression"}, false);
        if (type == "constant") value.solid.velocity_ic = new IC::Constant(value.geom, pp, "solid.velocity.ic.constant");
        else if (type == "expression") value.solid.velocity_ic = new IC::Expression(value.geom, pp, "solid.velocity.ic.expression");
    }
    {
        std::string type = "constant";
        // initial condition for momentum in the solid phase
        pp_query_validate("solid.momentum.ic.type", type, {"constant","expression"});
        if (type == "constant") value.solid.momentum_ic = new IC::Constant(value.geom, pp, "solid.momentum.ic.constant");
        else if (type == "expression") value.solid.momentum_ic = new IC::Expression(value.geom, pp, "solid.momentum.ic.expression");
    }
    {
        std::string type = "constant";
        // initial condition for energy in the solid phase
        pp_query_validate("solid.energy.ic.type", type, {"constant","expression"});
        if (type == "constant") value.solid.energy_ic = new IC::Constant(value.geom, pp, "solid.energy.ic.constant");
        else if (type == "expression") value.solid.energy_ic = new IC::Expression(value.geom, pp, "solid.energy.ic.expression");
    }
    {
        std::string type;
        pp_forbid("SolidDensity.ic.type","replaced by solid.density.ic.type");
        pp_forbid("SolidDensity.ic.constant","replaced by solid.density.ic.constant");
        pp_forbid("SolidDensity.ic.expression","replaced by solid.density.ic.expression");
        // initial condition for density in the solid phase
        pp_query_validate("solid.density.ic.type", type,{"constant","expression"});
        if (type == "constant") value.solid.density_ic = new IC::Constant(value.geom, pp, "solid.density.ic.constant");
        else if (type == "expression") value.solid.density_ic = new IC::Expression(value.geom, pp, "solid.density.ic.expression");
    }
    {
        std::string type;
        // initial condition for pressure in the solid phase
        pp_query_validate("solid.pressure.ic.type", type,{"constant","expression"});
        if (type == "constant") value.solid.pressure_ic = new IC::Constant(value.geom, pp, "solid.pressure.ic.constant");
        else if (type == "expression") value.solid.pressure_ic = new IC::Expression(value.geom, pp, "solid.pressure.ic.expression");
    }
    {
        std::string type;
        pp_forbid("Density.ic.type","replaced with density.ic.type");
        pp_forbid("Density.ic.constant","replaced with density.ic.constant");
        pp_forbid("Density.ic.expression","replaced with density.ic.expression");
        // initial condition for density
        pp_query_validate("density.ic.type", type, {"constant","expression"});
        if (type == "constant") value.density_ic = new IC::Constant(value.geom, pp, "density.ic.constant");
        else if (type == "expression") value.density_ic = new IC::Expression(value.geom, pp, "density.ic.expression");
    }
    {
        std::string type;
        // initial condition for density in interface region
        pp_query_validate("rhoInterface.ic.type", type,{"constant","expression"});
        if (type == "constant") value.ic_rhoInterface = new IC::Constant(value.geom, pp, "rhoInterface.ic.constant");
        else if (type == "expression") value.ic_rhoInterface = new IC::Expression(value.geom, pp, "rhoInterface.ic.expression");
        else Util::Abort(INFO, "Invalid rhoInterface.ic: ", type);
    }
    {
        std::string type = "constant";
        // injectied velocity initial condition
        pp_query_validate("vInjected.ic.type", type, {"constant","expression"});
        if (type == "constant") value.ic_vInjected = new IC::Constant(value.geom, pp, "vInjected.ic.constant");
        else if (type == "expression") value.ic_vInjected = new IC::Expression(value.geom, pp, "vInjected.ic.expression");
    }
    {
        std::string type;
        // heat flux initial condition
        pp_query_validate("q.ic.type", type,{"constant","expression"}); 
        if (type == "constant") value.ic_q = new IC::Constant(value.geom, pp, "q.ic.constant");
        else if (type == "expression") value.ic_q = new IC::Expression(value.geom, pp, "q.ic.expression");
    }

    if (pp.AnyUnusedInputs())
    {
        Util::Exception(INFO,"Exiting due to unused inputs");
    }
}


void Hydro::Initialize(int lev)
{
    BL_PROFILE("Integrator::Hydro::Initialize");
 
    eta_ic           ->Initialize(lev, eta_mf,     0.0);
    eta_ic           ->Initialize(lev, eta_old_mf, 0.0);
    etadot_mf[lev]   ->setVal(0.0);

    velocity_ic      ->Initialize(lev, velocity_mf, 0.0);
    pressure_ic      ->Initialize(lev, pressure_mf, 0.0);
    density_ic       ->Initialize(lev, density_mf, 0.0);

    solid.velocity_ic ->Initialize(lev, solid.velocity_mf, 0.0);
    solid.density_ic  ->Initialize(lev, solid.density_mf, 0.0);
    solid.pressure_ic ->Initialize(lev, solid.pressure_mf, 0.0);
    solid.momentum_ic ->Initialize(lev, solid.momentum_mf, 0.0);
    solid.energy_ic   ->Initialize(lev, solid.energy_mf, 0.0);

    ic_rhoInterface  ->Initialize(lev, rhoInterface_mf, 0.0);
    ic_vInjected     ->Initialize(lev, vInjected_mf,    0.0);
    ic_q             ->Initialize(lev, q_mf,            0.0);

    source_mf[lev]   ->setVal(0.0);

    Mix(lev);
}

void Hydro::Mix(int lev)
{
    Util::Message(INFO, eta_mf[lev]->nComp());
    for (amrex::MFIter mfi(*eta_mf[lev], false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        amrex::Array4<Set::Scalar> const& E_old = (*energy_old_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& E = (*energy_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& rho_old = (*density_old_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& rho = (*density_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& M_old = (*momentum_old_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& M = (*momentum_mf[lev]).array(mfi);

        amrex::Array4<const Set::Scalar> const& v = (*velocity_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& p = (*pressure_mf[lev]).array(mfi);

        amrex::Array4<const Set::Scalar> const& rho_solid = (*solid.density_mf[lev]).array(mfi);
        //amrex::Array4<const Set::Scalar> const& v_solid   = (*solid.velocity_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& M_solid   = (*solid.momentum_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& E_solid   = (*solid.energy_mf[lev]).array(mfi);

        amrex::Array4<const Set::Scalar> const& eta = (*eta_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {  
            rho(i, j, k) = eta(i, j, k) * rho(i, j, k) + (1.0 - eta(i, j, k)) * rho_solid(i, j, k);
            rho_old(i, j, k) = rho(i, j, k);

            E(i, j, k) =
                (
                    0.5 * (v(i, j, k, 0) * v(i, j, k, 0) + v(i, j, k, 1) * v(i, j, k, 1)) * rho(i, j, k)
                    +
                    (p(i, j, k) / (gamma - 1.0))
                    )
                * eta(i, j, k) 
                + 
                ( E_solid(i,j,k) ) * (1.0 - eta(i, j, k))
                ;
            E_old(i, j, k) = E(i, j, k);

            M(i, j, k, 0) = rho(i, j, k) * v(i, j, k, 0) * eta(i, j, k) + M_solid(i, j, k, 0) * (1.0 - eta(i, j, k));
            M(i, j, k, 1) = rho(i, j, k) * v(i, j, k, 1) * eta(i, j, k) + M_solid(i, j, k, 1) * (1.0 - eta(i, j, k));
            M_old(i, j, k, 0) = M(i, j, k, 0);
            M_old(i, j, k, 1) = M(i, j, k, 1);
        });
    }

    c_max = 0.0;
    vx_max = 0.0;
    vy_max = 0.0;
}

void Hydro::UpdateEta(int lev, Set::Scalar time)
{
    eta_ic->Initialize(lev, eta_mf, time);
}

void Hydro::TimeStepBegin(Set::Scalar, int /*iter*/)
{
}

void Hydro::TimeStepComplete(Set::Scalar, int lev)
{
    return;

    const Set::Scalar* DX = geom[lev].CellSize();

    amrex::ParallelDescriptor::ReduceRealMax(c_max);
    amrex::ParallelDescriptor::ReduceRealMax(vx_max);
    amrex::ParallelDescriptor::ReduceRealMax(vy_max);

    Set::Scalar new_timestep = cfl / ((c_max + vx_max) / DX[0] + (c_max + vy_max) / DX[1]);

    Util::Assert(INFO, TEST(AMREX_SPACEDIM == 2));

    SetTimestep(new_timestep);
}

#if AMREX_SPACEDIM==2
void Hydro::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    std::swap(eta_old_mf, eta_mf);
    std::swap(momentum_old_mf[lev], momentum_mf[lev]);
    std::swap(velocity_old_mf[lev], velocity_mf[lev]);
    std::swap(energy_old_mf[lev],   energy_mf[lev]);
    std::swap(density_old_mf[lev],  density_mf[lev]);
    std::swap(pressure_old_mf[lev], pressure_mf[lev]);

    UpdateEta(lev, time);
    const Set::Scalar* DX = geom[lev].CellSize();

    for (amrex::MFIter mfi(*eta_mf[lev], false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<const Set::Scalar> const& eta_new = (*eta_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& eta = (*eta_old_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar>       const& etadot = (*etadot_mf[lev]).array(mfi);

        Set::Patch<Set::Scalar>       omega = vorticity_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> v     = velocity_mf.Patch(lev,mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {

            etadot(i, j, k) = (eta_new(i, j, k) - eta(i, j, k)) / dt;

            auto sten = Numeric::GetStencil(i, j, k, bx);
            Set::Vector grad_ux = Numeric::Gradient(v, i, j, k, 0, DX, sten);
            Set::Vector grad_uy = Numeric::Gradient(v, i, j, k, 1, DX, sten);
            omega(i, j, k) = eta(i, j, k) * (grad_uy(0) - grad_ux(1));
        });
    }
    
    //
    // calculate source terms
    //
    for (amrex::MFIter mfi(*eta_mf[lev], false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        Set::Patch<const Set::Scalar>  rho = density_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar>  rho_solid = solid.density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar>  p_solid = solid.pressure_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar>  E   = energy_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> _M   = momentum_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> _u   = velocity_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar>  p   = pressure_old_mf.Patch(lev,mfi);

        Set::Patch<const Set::Scalar> eta    = eta_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> etadot = etadot_mf.Patch(lev,mfi);

        Set::Patch<const Set::Scalar> _M_solid = solid.momentum_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar>       source = source_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            //
            // Calculate gradients
            //
            Set::Vector u            = Set::Vector(_u(i,j,k,0),_u(i,j,k,1));
            Set::Matrix grad_u       = Numeric::Gradient(_u, i, j, k, DX);
            Set::Scalar div_u        = Numeric::Divergence(_u, i, j, k, DX);

            Set::Vector M            = Set::Vector(_M(i,j,k,0),_M(i,j,k,1));
            Set::Vector M_solid      = Set::Vector(_M_solid(i,j,k,0),_M_solid(i,j,k,1));
            
            Set::Vector grad_eta     = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();
            Set::Matrix hess_eta     = Numeric::Hessian(eta, i, j, k, 0, DX);
            
            //
            // Mass conservation source term
            //
            Set::Scalar mdot0 = 0.0; 
            // mixture terms
            mdot0 += -rho_solid(i,j,k) * grad_eta.dot(u);
            mdot0 += (1.0 - eta(i,j,k)) * rho_solid(i,j,k) * div_u;
            

            //
            // Momentum conservation source term
            //
            Set::Vector Pdot0 = Set::Vector::Zero();
            // passive terms
            Pdot0 += ( rho(i,j,k)*(u*u.transpose()) + mu*grad_u + 0.5*mu*grad_u.transpose() + 0.5*mu*div_u*Set::Matrix::Identity()) * grad_eta;
            // active terms [ none ftm ]
            // mixture terms - rho
            Pdot0 += - rho_solid(i,j,k) * (u * u.transpose()) * grad_eta; 
            Pdot0 += (1.0-eta(i,j,k)) * rho_solid(i,j,k) * grad_u * u;
            Pdot0 += (1.0-eta(i,j,k)) * rho_solid(i,j,k) * u * div_u;
            // mixture terms - pressure
            Pdot0 += - grad_eta * p_solid(i,j,k);

            //Set::Scalar etadot_cell  = etadot(i, j, k);
            //Set::Matrix gradu        = Numeric::Gradient(v, i, j, k, DX);
            // Set::Scalar divu         = 0.0;//Numeric::Divergence(v, i, j, k, 2, DX);
            //Set::Matrix I            = Set::Matrix::Identity();
        
            // Flow values
            //Set::Vector u(v(i,j,k,0),v(i,j,k,1));
            // Prescribed values
            //Set::Scalar rho0  = rhoInterface(i, j, k);
            //Set::Vector u0    = Set::Vector(vInjected(i,j,k,0),vInjected(i,j,k,1)) - grad_eta * etadot(i, j, k)/(grad_eta_mag * grad_eta_mag + 1.0e-12);
            Set::Vector u0 = Set::Vector::Zero(); // hard-coding imposed velocity for now
            //Set::Vector q0    = Set::Vector(q(i,j,k,0),q(i,j,k,1));
            //Set::Matrix T     = 0.5 * mu * (gradu + gradu.transpose());

            // Set::Scalar mdot0 = (rho0 * u0).dot(grad_eta);
            // Set::Vector Pdot0 = Set::Vector::Zero();//(rho0 * (u0*u0.transpose()) - T + 0.5 * (gradu.transpose() + divu * I))*grad_eta;
            // Set::Scalar qdot0 = (/*0.5*rho0*(u0.dot(u0))*u0 - T*u0 +*/ q0).dot(grad_eta); 
            Set::Vector Ldot0 = Set::Vector::Zero();

            for (int p = 0; p<2; p++)
            for (int q = 0; q<2; q++)
            for (int r = 0; r<2; r++)
            for (int s = 0; s<2; s++)
            {
                Set::Scalar Mpqrs = 0.0;
                if (p==r && q==s) Mpqrs += 0.5 * mu;
                if (p==s && q==r) Mpqrs += 0.5 * mu;
                Ldot0(p) += Mpqrs * (u(r) - u0(r)) * hess_eta(q, s);
            }
            
            Set::Scalar qdot0 = 0.0;

            source(i,j, k, 0) = mdot0;
            source(i,j, k, 1) = (Pdot0(0) + Ldot0(0));
            source(i,j, k, 2) = (Pdot0(1) + Ldot0(1));
            source(i,j, k, 3) = (qdot0);

            //Source Term Update
            // E(i, j, k)    += source(i, j, k, 3) * dt;
            // rho(i, j, k)  += source(i, j, k, 0) * dt;
            // M(i, j, k, 0) += source(i, j, k, 1) * dt;
            // M(i, j, k, 1) += source(i, j, k, 2) * dt;

            
        });
    }
    //source_mf[lev] -> FillBoundary();
    //return;
    

    for (amrex::MFIter mfi(*eta_mf[lev], false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        Set::Patch<const Set::Scalar> E   = energy_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> rho = density_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M   = momentum_old_mf.Patch(lev,mfi);

        Set::Patch<const Set::Scalar> eta    = eta_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> eta_new= eta_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> etadot = etadot_mf.Patch(lev,mfi);

        Set::Patch<const Set::Scalar> v = velocity_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> p = pressure_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> p_solid = solid.pressure_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar> E_new   = energy_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> rho_new = density_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> M_new   = momentum_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> v_new   = velocity_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> p_new   = pressure_mf.Patch(lev,mfi);


        Set::Patch<const Set::Scalar> source   = source_mf.Patch(lev,mfi);


        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            //Viscous Terms
            Set::Scalar lap_ux = Numeric::Laplacian(v, i, j, k, 0, DX);
            Set::Scalar lap_uy = Numeric::Laplacian(v, i, j, k, 1, DX);

            Set::Vector grad_eta     = Numeric::Gradient(eta, i, j, k, 0, DX);
            //Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();

            //Godunov flux
            Solver::Local::Riemann::Roe::State state_x(rho(i, j, k), v(i, j, k, 0), v(i, j, k, 1), p(i, j, k));
            Solver::Local::Riemann::Roe::State state_y(rho(i, j, k), v(i, j, k, 1), v(i, j, k, 0), p(i, j, k));

            Solver::Local::Riemann::Roe::State lo_statex(rho(i - 1, j, k), v(i - 1, j, k, 0), v(i - 1, j, k, 1), p(i - 1, j, k));
            Solver::Local::Riemann::Roe::State hi_statex(rho(i + 1, j, k), v(i + 1, j, k, 0), v(i + 1, j, k, 1), p(i + 1, j, k));

            Solver::Local::Riemann::Roe::State lo_statey(rho(i, j - 1, k), v(i, j - 1, k, 1), v(i, j - 1, k, 0), p(i, j - 1, k));
            Solver::Local::Riemann::Roe::State hi_statey(rho(i, j + 1, k), v(i, j + 1, k, 1), v(i, j + 1, k, 0), p(i, j + 1, k));

            Solver::Local::Riemann::Roe::Flux flux_xlo, flux_ylo, flux_xhi, flux_yhi;

            try
            {
                //lo interface fluxes
                flux_xlo = Solver::Local::Riemann::Roe::Solve(lo_statex, state_x, gamma);
                flux_ylo = Solver::Local::Riemann::Roe::Solve(lo_statey, state_y, gamma);

                //hi interface fluxes
                flux_xhi = Solver::Local::Riemann::Roe::Solve(state_x, hi_statex, gamma);
                flux_yhi = Solver::Local::Riemann::Roe::Solve(state_y, hi_statey, gamma);
            }
            catch (...)
            {
                Util::Message(INFO,"i=",i," j=",j);
                Util::Exception(INFO,"riemann solver failed");
            }

            // Set::Scalar E_solid = (0.5 * (v_solid(i, j, k, 0) * v_solid(i, j, k, 0) + v_solid(i, j, k, 1) * v_solid(i, j, k, 1)) * rho_solid(i, j, k) + p(i, j, k) / (gamma - 1.0));


            rho_new(i, j, k) =
                rho(i,j,k)
                + 
                (
                    // advection terms
                    (flux_xlo.mass - flux_xhi.mass) / DX[0] +
                    (flux_ylo.mass - flux_yhi.mass) / DX[1]
                    //+
                    // source terms
                    + source(i,j,k,0)
                ) * dt 

                // Transient
                // - (rho(i, j, k) - rho_solid(i, j, k)) * etadot(i, j, k) * dt

                /*Update solid density*/
                //+ (1.0 - eta(i, j, k)) * (rho_solid(i, j, k) - rho(i,j,k))
                ;         
                
            M_new(i, j, k, 0) =
                M(i, j, k, 0)
                + //eta(i, j, k) * 
                (
                    // Advection terms
                    (flux_xlo.momentum_normal  - flux_xhi.momentum_normal) / DX[0] +
                    (flux_ylo.momentum_tangent - flux_yhi.momentum_tangent) / DX[1] 
                    +
                    // Viscous terms
                    (mu * lap_ux)
                    + source(i,j,k,1)
                )*dt
                // -  (M(i, j, k, 0) - rho_solid(i, j, k) * v_solid(i, j, k, 0)) * etadot(i, j, k) * dt

                /*Update solid momentum*/   
                //+ (1.0 - eta(i, j, k)) * (rho_solid(i, j, k) * v_solid(i, j, k, 0) - M(i,j,k,0) )
                ;
                
            M_new(i, j, k, 1) =
                M(i, j, k, 1)
                + //eta(i, j, k) * 
                (
                    // Advection terms
                    (flux_xlo.momentum_tangent - flux_xhi.momentum_tangent) / DX[0] +
                    (flux_ylo.momentum_normal  - flux_yhi.momentum_normal ) / DX[1]
                    +
                    // Viscous terms
                    (mu * lap_uy)
                    + source(i,j,k,2)
                ) * dt

                // moving eta
                //- (M(i, j, k, 1) - rho_solid(i, j, k) * v_solid(i, j, k, 1)) * etadot(i, j, k) * dt

                //*Update solid momentum*/   
                //+ (1.0 - eta(i, j, k)) * (rho_solid(i, j, k) * v_solid(i, j, k, 1) - M(i,j,k,1))
                ;

            //Godunov fluxes
            E_new(i, j, k) =
                E(i, j, k)
                + 
                /*Update fluid energy*/
                // eta(i, j, k) * 
                (
                    (flux_xlo.energy - flux_xhi.energy) / DX[0]
                    +
                    (flux_ylo.energy - flux_yhi.energy) / DX[1]
                    ///////////+ 2. * mu * (div_u * div_u + div_u * symgrad_u) - 2./3. * mu * div_u * div_u;
                    //+
                    //source(i,j,k,3)
                ) * dt 

                //Should be jump value
                // - (E(i, j, k) - E_solid) * etadot(i, j, k) * dt 

                /*Update solid energnsy*/
                //+ (1.0 - eta(i, j, k)) * (E_solid - E(i,j,k))
                ;


            v_new(i, j, k, 0) = eta_new(i,j,k) * M_new(i, j, k, 0) / rho_new(i, j, k);
            v_new(i, j, k, 1) = eta_new(i,j,k) * M_new(i, j, k, 1) / rho_new(i, j, k);
            p_new(i, j, k) =
                eta_new(i,j,k) * (
                    (E_new(i, j, k) - 0.5 * (M_new(i, j, k, 0) * M_new(i, j, k, 0) + M_new(i, j, k, 1) * M_new(i, j, k, 1))
                     / rho_new(i, j, k)) * (gamma - 1.0)
                    )
                +
                (1.0 - eta_new(i,j,k)) * p_solid(i,j,k);

        });
    }
}//end Advance
#else
void Hydro::Advance(int, Set::Scalar, Set::Scalar)
{
    Util::Abort(INFO,"Works in 2D only");
}
#endif




void Hydro::Regrid(int lev, Set::Scalar /* time */)
{
    BL_PROFILE("Integrator::Hydro::Regrid");
    if (lev < finest_level) return;

    Util::Message(INFO, "Regridding on level", lev);
}//end regrid

//void Hydro::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar time, int ngrow)
void Hydro::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar, int)
{
    BL_PROFILE("Integrator::Flame::TagCellsForRefinement");

    const Set::Scalar* DX = geom[lev].CellSize();
    Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

    // Eta criterion for refinement
    for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<const Set::Scalar> const& eta = (*eta_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
            if (grad_eta.lpNorm<2>() * dr * 2 > eta_refinement_criterion) tags(i, j, k) = amrex::TagBox::SET;
        });
    }

    // Vorticity criterion for refinement
    for (amrex::MFIter mfi(*vorticity_mf[lev], true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<const Set::Scalar> const& omega = (*vorticity_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Vector grad_omega = Numeric::Gradient(omega, i, j, k, 0, DX);
            if (grad_omega.lpNorm<2>() * dr * 2 > omega_refinement_criterion) tags(i, j, k) = amrex::TagBox::SET;
        });
    }

}//end TagCells

// void Hydro::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/, const amrex::MFIter &mfi, const amrex::Box &box)
// {
//   BL_PROFILE("Hydro::Integrate");
//   const Set::Scalar *DX = geom[amrlev].CellSize();
//   Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
//   amrex::Array4<amrex::Real> const &eta = (*eta_mf[amrlev]).array(mfi);
//   amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k){
//   volume += eta(i, j, k, 0) * dv;
//   Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
//   Set::Scalar normgrad = grad.lpNorm<2>();
//   Set::Scalar da = normgrad * dv;
//   area += da;
//   });


//  }//end Integrate

//  void Hydro::UpdateModel(int /*a_step*/)
//  {
//    for (int lev = 0; lev <= finest_level; ++lev)
//    {
//      eta_mf[lev] -> FillBoundary();
//      density_mf[lev] -> FillBoundary();
//      energy_mf[lev] -> FillBoundary();
//      Momentum[lev] -> FillBoundary();
//      
//      for (MFIter mfi(*model_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
//      {
//  amrex::Box bx = mfi.nodaltilebox();
//  amrex::Array4<model_type> const &model = model_mf[lev]->array(mfi);
//
//  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
//  // TODO
//
//  });
//
//
//     } // end For2

//     Util::RealFillBoundary(*model_mf[lev], geom[lev]);
//     amrex::MultiFab::Copy(*psi_mf[lev], *eta_mf[lev], 0, 0, 1, psi_mf[lev]-> nGrow());

//    } //end For1
//}//end update


}//end code
