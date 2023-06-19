/*
 * Cell_SPM_diffusion.cpp
 *
 * Implements the functions for the parent class of the Cells
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "Cell_SPM.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <array>
#include <algorithm>
#include <utility>

namespace slide {
ValuePair<double> Cell_SPM::calcSurfaceConcentration(ValuePair<double> j, ValuePair<double> Dt) //!< Should not throw normally, except divide by zero?
{
  using settings::nch;
  ValuePair<double> c_surf{ 0, 0 }; //!< positive and negative surface concentrations.
  //!< Calculate the surface concentration at the positive particle
  //!< 	c_surf.pos = M->Cp[0][:] * zp[:] + M->Dp*jp/Dpt
  for (size_t i = 0; i < nch; i++)
    c_surf.pos += M->Cp[0][i] * st.zp(i);
  c_surf.pos += M->Dp[0] * j.pos / Dt.pos;

  //!< Calculate the surface concentration at the negative particle
  //!< 	c_surf.neg = M->Cn[0][:] * zn[:] + M->Dn*jn/Dnt
  for (size_t i = 0; i < nch; i++)
    c_surf.neg += M->Cn[0][i] * st.zn(i);
  c_surf.neg += M->Dn[0] * j.neg / Dt.neg;

  return c_surf;
}

ValuePair<> Cell_SPM::calcOverPotential(ValuePair<> cs, double i_app) //!< Should not throw normally, except divide by zero?
{
  using namespace PhyConst;

  const auto ArrheniusCoeff = calcArrheniusCoeff();

  //!< Calculate the rate constants at the cell's temperature using an Arrhenius relation
  const double kpt = kp * std::exp(kp_T * ArrheniusCoeff); //!< Rate constant at the positive electrode at the cell's temperature [m s-1]
  const double knt = kn * std::exp(kn_T * ArrheniusCoeff); //!< Rate constant at the negative electrode at the cell's temperature [m s-1]

  //!< Calculate the overpotential using the Bulter-Volmer equation
  //!< if alpha is 0.5, the Bulter-Volmer relation can be inverted to eta = 2RT / (nF) asinh(x)
  //!< and asinh(x) = ln(x + sqrt(1+x^2) -> to asinh(x) function.
  const double i0p = kpt * n * F * std::sqrt(C_elec * cs.pos * (Cmaxpos - cs.pos)); //!< exchange current density of the positive electrode
  const double i0n = knt * n * F * std::sqrt(C_elec * cs.neg * (Cmaxneg - cs.neg)); //!< exchange current density of the negative electrode
  const double xp = -0.5 * i_app / (st.ap() * st.thickp()) / i0p;                   //!< x for the cathode
  const double xn = 0.5 * i_app / (st.an() * st.thickn()) / i0n;                    //!< x for the anode
  const double etap = (2 * Rg * st.T()) / (n * F) * std::asinh(xp);                 //!< cathode overpotential [V], < 0 on discharge
  const double etan = (2 * Rg * st.T()) / (n * F) * std::asinh(xn);                 //!< anode overpotential [V],  > 0 on discharge

  return { etap, etan };
}

std::array<double, 3> Cell_SPM::calcMolarFlux() //!< Should not throw normally, except divide by zero?
{
  using PhyConst::F;
  const double i_app = I() / geo.elec_surf;                   //!< current density on the electrode [A m-2]
  const double jp = -i_app / (st.ap() * n * F * st.thickp()); //!< molar flux on the positive particle [mol m-2 s-1]
  const double jn = i_app / (st.an() * n * F * st.thickn());  //!< molar flux on the negative particle [mol m-2 s-1]

  return { i_app, jp, jn };
}

ValuePair<double> Cell_SPM::calcDiffConstant() //!< Should not throw normally, except divide by zero?
{
  const auto ArrheniusCoeff = calcArrheniusCoeff();
  //!< Calculate the diffusion constant at the battery temperature using an Arrhenius relation
  const double Dpt = st.Dp() * std::exp(Dp_T * ArrheniusCoeff); //!< diffusion constant of the positive particle [m s-1]
  const double Dnt = st.Dn() * std::exp(Dn_T * ArrheniusCoeff); //!< diffusion constant of the negative particle [m s-1]

  return { Dpt, Dnt };
}
} // namespace slide