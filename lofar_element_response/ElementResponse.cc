//# ElementResponse.cc: Functions to compute the (idealized) response of a LOFAR
//# LBA or HBA dual dipole antenna.
//#
//# Copyright (C) 2011
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$

#include "ElementResponse.h"
#include <cmath>

// The coefficients are kept in an unnamed namespace which effectively makes
// them invisible outside this translation unit.
namespace
{
#include "DefaultCoeffLBA.cc"
#include "DefaultCoeffHBA.cc"
}

namespace LOFAR
{

void element_response_lba(double freq, double az, double el,
    std::complex<double> (&response)[2][2])
{
    element_response(freq, az, el, response, default_lba_freq_center,
        default_lba_freq_range, default_lba_coeff_shape, default_lba_coeff);
}

void element_response_hba(double freq, double az, double el,
    std::complex<double> (&response)[2][2])
{
    element_response(freq, az, el, response, default_hba_freq_center,
        default_hba_freq_range, default_hba_coeff_shape, default_hba_coeff);
}

void element_response(double freq, double az, double el,
    std::complex<double> (&response)[2][2], double freq_center,
    double freq_range, const unsigned int (&coeff_shape)[3],
    const std::complex<double> coeff[])
{
    const unsigned int nHarmonics  = coeff_shape[0];
    const unsigned int nPowerTheta = coeff_shape[1];
    const unsigned int nPowerFreq  = coeff_shape[2];

    // The model is parameterized in terms of zenith angle. The appropriate
    // conversion is taken care of below.
    const double theta = (M_PI / 2) - el;

    // The model calculates azimuth from the +x dipole (SW direction)
    // LOFAR uses azimuth eastward positive from north.
    const double phi = az - (5 * M_PI / 4);

    // The model is parameterized in terms of a normalized frequency in the
    // range [-1, 1]. The appropriate conversion is taken care of below.
    freq = (freq - freq_center) / freq_range;

    // Initialize the response to zero.
    response[0][0] = 0.0;
    response[0][1] = 0.0;
    response[1][0] = 0.0;
    response[1][1] = 0.0;

    // The variables sign and kappa are used to compute the value of kappa
    // mentioned in the description of the beam model [kappa = (-1)^k * (2 * k
    //+ 1)] incrementally.
    int sign = 1, kappa = 1;

    std::complex<double> P[2], Pj[2];
    for(unsigned int k = 0; k < nHarmonics; ++k)
    {
        // Compute the (diagonal) projection matrix P for the current harmonic.
        // This requires the evaluation of two polynomials in theta and freq (of
        // degree nPowerTheta in theta and nPowerFreq in freq), one for each
        // element of P. The polynomials are evaluated using Horner's rule.

        // Horner's rule requires backward iteration of the coefficients, so
        // move the iterator to the first location past the end of the block of
        // coefficients for the current harmonic (k).
        coeff += nPowerTheta * nPowerFreq * 2;

        // Evaluate the polynomial. Note that the iterator is always decremented
        // before it is dereferenced, using the prefix operator--. After
        // evaluation of the polynomial, the iterator points exactly to the
        // beginning of the block of coefficients for the current harmonic (k),
        // that is, all the decrements together exactly cancel the increment
        // aplied above.

        // Evaluate the highest order term. Note that the order of the
        // assigments is important because of the iterator decrement, i.e. P[1]
        // should be assigned first.
        P[1] = *--coeff;
        P[0] = *--coeff;

        for(unsigned int i = 0; i < nPowerFreq - 1; ++i)
        {
            P[1] = P[1] * freq + *--coeff;
            P[0] = P[0] * freq + *--coeff;
        }

        // Evaluate the remaining terms.
        for(unsigned int j = 0; j < nPowerTheta - 1; ++j)
        {
            Pj[1] = *--coeff;
            Pj[0] = *--coeff;

            for(unsigned int i = 0; i < nPowerFreq - 1; ++i)
            {
                Pj[1] = Pj[1] * freq + *--coeff;
                Pj[0] = Pj[0] * freq + *--coeff;
            }

            P[1] = P[1] * theta + Pj[1];
            P[0] = P[0] * theta + Pj[0];
        }

        // Because the increment and decrements cancel, the iterator points to
        // the same location as at the beginning of this iteration of the outer
        // loop. The next iteration should use the coefficients for the next
        // harmonic (k), so we move the iterator to the start of that block.
        coeff += nPowerTheta * nPowerFreq * 2;

        // Compute the Jones matrix for the current harmonic, by rotating P over
        // kappa * phi, and add it to the result.
        const double angle = sign * kappa * phi;
        const double cp = std::cos(angle);
        const double sp = std::sin(angle);

        response[0][0] += cp * P[0];
        response[0][1] += -sp * P[1];
        response[1][0] += sp * P[0];
        response[1][1] += cp * P[1];

        // Update sign and kappa.
        sign = -sign;
        kappa += 2;
    }
}

} //# namespace lofar
