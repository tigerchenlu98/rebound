/**
 * @file 	integrator.c
 * @brief 	BS integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements the Gragg-Bulirsch-Stoer integration scheme.  
 *          It is a reimplementation of the fortran code by E. Hairer and G. Wanner.
 *          The starting point was the JAVA implementation in hipparchus:
 *          https://github.com/Hipparchus-Math/hipparchus/blob/master/hipparchus-ode/src/main/java/org/hipparchus/ode/nonstiff/GraggBulirschStoerIntegrator.java
 *
 * @section 	LICENSE
 * Copyright (c) 2021 Hanno Rein
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright (c) 2004, Ernst Hairer
 *
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the following
 * conditions are met:
 * 
 *  - Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 * BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <float.h> // for DBL_MAX
#include "rebound.h"
#include "integrator_bs.h"
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

void setStabilityCheck(struct reb_simulation_integrator_bs* ri_bs, int performStabilityCheck, int maxNumIter, int maxNumChecks, double stepsizeReductionFactor) {
    ri_bs->performTest = performStabilityCheck;
    ri_bs->maxIter     = (maxNumIter   <= 0) ? 2 : maxNumIter;
    ri_bs->maxChecks   = (maxNumChecks <= 0) ? 1 : maxNumChecks;

    if ((stepsizeReductionFactor < 0.0001) || (stepsizeReductionFactor > 0.9999)) {
        ri_bs->stabilityReduction = 0.5;
    } else {
        ri_bs->stabilityReduction = stepsizeReductionFactor;
    }
}

void setControlFactors(struct reb_simulation_integrator_bs* ri_bs, double control1, double control2, double control3, double control4) {

    if ((control1 < 0.0001) || (control1 > 0.9999)) {
        ri_bs->stepControl1 = 0.65;
    } else {
        ri_bs->stepControl1 = control1;
    }

    if ((control2 < 0.0001) || (control2 > 0.9999)) {
        ri_bs->stepControl2 = 0.94;
    } else {
        ri_bs->stepControl2 = control2;
    }

    if ((control3 < 0.0001) || (control3 > 0.9999)) {
        ri_bs->stepControl3 = 0.02;
    } else {
        ri_bs->stepControl3 = control3;
    }

    if ((control4 < 1.0001) || (control4 > 999.9)) {
        ri_bs->stepControl4 = 4.0;
    } else {
        ri_bs->stepControl4 = control4;
    }

}

void initializeArrays(struct reb_simulation_integrator_bs* ri_bs) {

    int size = ri_bs->maxOrder / 2;

    if ((ri_bs->sequence == NULL) || (ri_bs->sequence_length != size)) {
        // all arrays should be reallocated with the right size
        ri_bs->sequence        = realloc(ri_bs->sequence,sizeof(int)*size);
        ri_bs->costPerStep     = realloc(ri_bs->costPerStep,sizeof(int)*size);
        ri_bs->coeff           = realloc(ri_bs->coeff,sizeof(double*)*size);
        for (int k = ri_bs->sequence_length; k < size; ++k) {
            ri_bs->coeff[k] = NULL;
        }
        ri_bs->costPerTimeUnit = realloc(ri_bs->costPerTimeUnit,sizeof(double)*size);
        ri_bs->optimalStep     = realloc(ri_bs->optimalStep,sizeof(double)*size);
    }

    // step size sequence: 2, 6, 10, 14, ...
    for (int k = 0; k < size; ++k) {
        ri_bs->sequence[k] = 4 * k + 2;
    }

    // initialize the order selection cost array
    // (number of function calls for each column of the extrapolation table)
    ri_bs->costPerStep[0] = ri_bs->sequence[0] + 1;
    for (int k = 1; k < size; ++k) {
        ri_bs->costPerStep[k] = ri_bs->costPerStep[k - 1] + ri_bs->sequence[k];
    }

    // initialize the extrapolation tables
    for (int k = 1; k < size; ++k) {
        ri_bs->coeff[k] = realloc(ri_bs->coeff[k],sizeof(double)*k);
        for (int l = 0; l < k; ++l) {
            double ratio = ((double) ri_bs->sequence[k]) / ri_bs->sequence[k - l - 1];
            ri_bs->coeff[k][l] = 1.0 / (ratio * ratio - 1.0);
        }
    }

}

void setInterpolationControl(struct reb_simulation_integrator_bs* ri_bs, int useInterpolationErrorForControl, int mudifControlParameter) {

    ri_bs->useInterpolationError = useInterpolationErrorForControl;

    if ((mudifControlParameter <= 0) || (mudifControlParameter >= 7)) {
        ri_bs->mudif = 4;
    } else {
        ri_bs->mudif = mudifControlParameter;
    }

}

void setOrderControl(struct reb_simulation_integrator_bs* ri_bs, int maximalOrder, double control1, double control2) {

    if (maximalOrder > 6 && maximalOrder % 2 == 0) {
        ri_bs->maxOrder = maximalOrder;
    } else {
        ri_bs->maxOrder = 18;
    }

    if ((control1 < 0.0001) || (control1 > 0.9999)) {
        ri_bs->orderControl1 = 0.8;
    } else {
        ri_bs->orderControl1 = control1;
    }

    if ((control2 < 0.0001) || (control2 > 0.9999)) {
        ri_bs->orderControl2 = 0.9;
    } else {
        ri_bs->orderControl2 = control2;
    }

    // reinitialize the arrays
    initializeArrays(ri_bs);

}

void bs_constructor(struct reb_simulation* r){
    struct reb_simulation_integrator_bs* ri_bs = &(r->ri_bs);
    ri_bs->minStep = fabs(ri_bs->minStep);
    ri_bs->maxStep = fabs(ri_bs->maxStep);
    ri_bs->initialStep = -1;
    setStabilityCheck(ri_bs, 1, -1, -1, -1);
    setControlFactors(ri_bs, -1, -1, -1, -1);
    setOrderControl(ri_bs, -1, -1, -1);
    setInterpolationControl(ri_bs, 1, -1);
}

double* computeDerivatives(double t, double* const yEnd){
    // DUMMY
    return 0; // TODO 
}


int tryStep(struct reb_simulation_integrator_bs* ri_bs, const double t0, const double* y0, const int y0_length, const double step, const int k, const double* scale, const int scale_length, double** const f, double* const yMiddle, double* const yEnd) {

    const int    n        = ri_bs->sequence[k];
    const double subStep  = step / n;
    const double subStep2 = 2 * subStep;

    // first substep
    double t = t0 + subStep;
    for (int i = 0; i < y0_length; ++i) {
        yEnd[i] = y0[i] + subStep * f[0][i];
    }
    f[1] = computeDerivatives(t, yEnd);

    // other substeps
    double* const yTmp = malloc(sizeof(double)*y0_length); // IMPROVE: should allocate this only once
    for (int j = 1; j < n; ++j) {

        if (2 * j == n) {
            // save the point at the middle of the step
            for (int i = 0; i < y0_length; ++i) {
                yMiddle[i] = yEnd[i];
            }
        }

        t += subStep;
        for (int i = 0; i < y0_length; ++i) {
            const double middle = yEnd[i];
            yEnd[i]       = y0[i] + subStep2 * f[j][i];
            yTmp[i]       = middle;
        }

        f[j + 1] = computeDerivatives(t, yEnd);

        // stability check
        if (ri_bs->performTest && (j <= ri_bs->maxChecks) && (k < ri_bs->maxIter)) {
            double initialNorm = 0.0;
            for (int l = 0; l < scale_length; ++l) {
                const double ratio = f[0][l] / scale[l];
                initialNorm += ratio * ratio;
            }
            double deltaNorm = 0.0;
            for (int l = 0; l < scale_length; ++l) {
                const double ratio = (f[j + 1][l] - f[0][l]) / scale[l];
                deltaNorm += ratio * ratio;
            }
            if (deltaNorm > 4 * MAX(1.0e-15, initialNorm)) {
                return 0;
            }
        }

    }

    // correction of the last substep (at t0 + step)
    for (int i = 0; i < y0_length; ++i) {
        yEnd[i] = 0.5 * (yTmp[i] + yEnd[i] + subStep * f[n][i]);
    }

    free(yTmp);
    return 1;

}

void extrapolate(struct reb_simulation_integrator_bs* ri_bs, const int offset, const int k, double** const diag, double* const last, const int last_length) {
    // update the diagonal
    for (int j = 1; j < k; ++j) {
        for (int i = 0; i < last_length; ++i) {
            // Aitken-Neville's recursive formula
            diag[k - j - 1][i] = diag[k - j][i] +
                ri_bs->coeff[k + offset][j - 1] * (diag[k - j][i] - diag[k - j - 1][i]);
        }
    }

    // update the last element
    for (int i = 0; i < last_length; ++i) {
        // Aitken-Neville's recursive formula
        last[i] = diag[0][i] + ri_bs->coeff[k + offset][k - 1] * (diag[0][i] - last[i]);
    }
}

double ulp(double x){
    return nextafter(x, INFINITY) - x;
}


struct ExpandableODE{

};
struct ODEState{
    double t; // getTime()
    double* y; // primary state
    int y_length;
    double* yDot; // secondary state
};

void sanityChecks(const struct ODEState initialState, const double t) {
    const double threshold = 1000 * ulp(MAX(fabs(initialState.t), fabs(t)));
    const double dt = fabs(initialState.t - t);
    if (dt <= threshold) {
        printf("Error. Integration interval too small.");
        exit(0);
    }
    // TODO set mainSetDimension

}

double getTolerance(struct reb_simulation_integrator_bs* ri_bs, int i, double scale){
    return ri_bs->scalAbsoluteTolerance + ri_bs->scalRelativeTolerance * scale;
}

void rescale(struct reb_simulation_integrator_bs* ri_bs, double* const y1, double* const y2, double* const scale, int scale_length) {
    for (int i = 0; i < scale_length; ++i) {
        scale[i] = getTolerance(ri_bs, i, MAX(fabs(y1[i]), fabs(y2[i])));
    }
} 
double filterStep(struct reb_simulation_integrator_bs* ri_bs, const double h, const int forward, const int acceptSmall){
    double filteredH = h;
    if (fabs(h) < ri_bs->minStep) {
        if (acceptSmall) {
            filteredH = forward ? ri_bs->minStep : -ri_bs->minStep;
        } else {
            printf("Error. Minimal stepsize reached during integration.");
            exit(0);
        }
    }

    if (filteredH > ri_bs->maxStep) {
        filteredH = ri_bs->maxStep;
    } else if (filteredH < -ri_bs->maxStep) {
        filteredH = -ri_bs->maxStep;
    }

    return filteredH;

}

struct ODEState integrate(struct reb_simulation* r, struct reb_simulation_integrator_bs* ri_bs, const struct ExpandableODE equations, const struct ODEState initialState, const double finalTime){

    sanityChecks(initialState, finalTime);
    // TODO setStepStart 
    const int forward = finalTime > initialState.t;

    // create some internal working arrays
    int y_length = initialState.y_length;
    double*        y         = initialState.y;
    double* const  y1        = malloc(sizeof(double)*y_length); // TODO free
    double** const diagonal  = malloc(sizeof(double*)*(ri_bs->sequence_length - 1)); // TODO free
    double** const y1Diag    = malloc(sizeof(double*)*(ri_bs->sequence_length - 1)); // TODO free
    for (int k = 0; k < ri_bs->sequence_length - 1; ++k) {
        diagonal[k] = malloc(sizeof(double)*y_length); // TODO free
        y1Diag[k]   = malloc(sizeof(double)*y_length); // TODO free
    }

    double*** const fk = malloc(sizeof(double**)*ri_bs->sequence_length); // TODO free
    for (int k = 0; k < ri_bs->sequence_length; ++k) {
        fk[k] = malloc(sizeof(double*)*(ri_bs->sequence[k] + 1)); // TODO free
    }

    // scaled derivatives at the middle of the step $\tau$
    // (element k is $h^{k} d^{k}y(\tau)/dt^{k}$ where h is step size...)
    double** const yMidDots = malloc(sizeof(double*)*(1 + 2 * ri_bs->sequence_length)); // TODO free
    for (int k = 0; k < ri_bs->sequence_length; ++k) {
        yMidDots[k] = malloc(sizeof(double*)*y_length); // TODO free
    }

    // initial scaling
    const int mainSetDimension = y_length; 
    double* const scale = malloc(sizeof(double)*mainSetDimension); // TODO free
    rescale(ri_bs, y, y, scale, y_length);

    // initial order selection
    const double tol    = ri_bs->scalRelativeTolerance;
    const double log10R = log10(MAX(1.0e-10, tol));
    int targetIter = MAX(1,
            MIN(ri_bs->sequence_length - 2,
                (int) floor(0.5 - 0.6 * log10R)));

    double  hNew                     = 0;
    double  maxError                 = DBL_MAX;
    int previousRejected         = 0;
    int firstTime                = 1;
    int newStep                  = 1;
    ri_bs->costPerTimeUnit[0] = 0;
    ri_bs->isLastStep = 0;
    do {

        double error;
        int reject = 0;

        if (newStep) {

            // first evaluation, at the beginning of the step
            double* const yDot0 = initialState.yDot;
            for (int k = 0; k < ri_bs->sequence_length; ++k) {
                // all sequences start from the same point, so we share the derivatives
                fk[k][0] = yDot0;
            }

            if (firstTime) {
                hNew = r->dt; // Was able to guess step size: initializeStep(forward, 2 * targetIter + 1, scale, getStepStart(), equations.getMapper());
            }

            newStep = 0;

        }

        ri_bs->stepSize = hNew;

        // step adjustment near bounds
        if (forward) {
            if (initialState.t + ri_bs->stepSize >= finalTime) {
                ri_bs->stepSize = finalTime - initialState.t;
            }
        } else {
            if (initialState.t + ri_bs->stepSize <= finalTime) {
                ri_bs->stepSize = finalTime - initialState.t;
            }
        }
        const double nextT = initialState.t + ri_bs->stepSize;
        ri_bs->isLastStep = (forward ? (nextT >= finalTime) : (nextT <= finalTime));

        // iterate over several substep sizes
        int k = -1;
        for (int loop = 1; loop; ) {

            ++k;

            // modified midpoint integration with the current substep
            if ( ! tryStep(ri_bs, initialState.t, y, y_length, ri_bs->stepSize, k, scale, y_length, fk[k],
                        (k == 0) ? yMidDots[0] : diagonal[k - 1],
                        (k == 0) ? y1 : y1Diag[k - 1])) {

                // the stability check failed, we reduce the global step
                hNew   = fabs(filterStep(ri_bs, ri_bs->stepSize * ri_bs->stabilityReduction, forward, 0));
                reject = 1;
                loop   = 0;

            } else {

                // the substep was computed successfully
                if (k > 0) {

                    // extrapolate the state at the end of the step
                    // using last iteration data
                    extrapolate(ri_bs, 0, k, y1Diag, y1, y_length);
                    rescale(ri_bs, y, y1, scale, y_length);

                    // estimate the error at the end of the step.
                    error = 0;
                    for (int j = 0; j < mainSetDimension; ++j) {
                        const double e = fabs(y1[j] - y1Diag[0][j]) / scale[j];
                        error += e * e;
                    }
                    error = sqrt(error / mainSetDimension);
                    if (isnan(error)) {
                        printf("Error. NaN appearing during integration.");
                        exit(0);
                    }

                    if ((error > 1.0e15) || ((k > 1) && (error > maxError))) {
                        // error is too big, we reduce the global step
                        hNew   = fabs(filterStep(ri_bs, ri_bs->stepSize * ri_bs->stabilityReduction, forward, 0));
                        reject = 1;
                        loop   = 0;
                    } else {

                        maxError = MAX(4 * error, 1.0);

                        // compute optimal stepsize for this order
                        const double exp = 1.0 / (2 * k + 1);
                        double fac = ri_bs->stepControl2 / pow(error / ri_bs->stepControl1, exp);
                        const double power = pow(ri_bs->stepControl3, exp);
                        fac = MAX(power / ri_bs->stepControl4, MIN(1. / power, fac));
                        const int acceptSmall = k < targetIter;
                        ri_bs->optimalStep[k]     = fabs(filterStep(ri_bs, ri_bs->stepSize * fac, forward, acceptSmall));
                        ri_bs->costPerTimeUnit[k] = ri_bs->costPerStep[k] / ri_bs->optimalStep[k];

                        // check convergence
                        switch (k - targetIter) {

                            case -1 :
                                if ((targetIter > 1) && ! previousRejected) {

                                    // check if we can stop iterations now
                                    if (error <= 1.0) {
                                        // convergence have been reached just before targetIter
                                        loop = 0;
                                    } else {
                                        // estimate if there is a chance convergence will
                                        // be reached on next iteration, using the
                                        // asymptotic evolution of error
                                        const double ratio = ((double) ri_bs->sequence[targetIter] * ri_bs->sequence[targetIter + 1]) / (ri_bs->sequence[0] * ri_bs->sequence[0]);
                                        if (error > ratio * ratio) {
                                            // we don't expect to converge on next iteration
                                            // we reject the step immediately and reduce order
                                            reject = 1;
                                            loop   = 0;
                                            targetIter = k;
                                            if ((targetIter > 1) &&
                                                    (ri_bs->costPerTimeUnit[targetIter - 1] <
                                                     ri_bs->orderControl1 * ri_bs->costPerTimeUnit[targetIter])) {
                                                --targetIter;
                                            }
                                            hNew = filterStep(ri_bs, ri_bs->optimalStep[targetIter], forward, 0);
                                        }
                                    }
                                }
                                break;

                            case 0:
                                if (error <= 1.0) {
                                    // convergence has been reached exactly at targetIter
                                    loop = 0;
                                } else {
                                    // estimate if there is a chance convergence will
                                    // be reached on next iteration, using the
                                    // asymptotic evolution of error
                                    const double ratio = ((double) ri_bs->sequence[k + 1]) / ri_bs->sequence[0];
                                    if (error > ratio * ratio) {
                                        // we don't expect to converge on next iteration
                                        // we reject the step immediately
                                        reject = 1;
                                        loop = 0;
                                        if ((targetIter > 1) &&
                                                (ri_bs->costPerTimeUnit[targetIter - 1] <
                                                 ri_bs->orderControl1 * ri_bs->costPerTimeUnit[targetIter])) {
                                            --targetIter;
                                        }
                                        hNew = filterStep(ri_bs, ri_bs->optimalStep[targetIter], forward, 0);
                                    }
                                }
                                break;

                            case 1 :
                                if (error > 1.0) {
                                    reject = 1;
                                    if ((targetIter > 1) &&
                                            (ri_bs->costPerTimeUnit[targetIter - 1] <
                                             ri_bs->orderControl1 * ri_bs->costPerTimeUnit[targetIter])) {
                                        --targetIter;
                                    }
                                    hNew = filterStep(ri_bs, ri_bs->optimalStep[targetIter], forward, 0);
                                }
                                loop = 0;
                                break;

                            default :
                                if ((firstTime || ri_bs->isLastStep) && (error <= 1.0)) {
                                    loop = 0;
                                }
                                break;

                        }

                    }
                }
            }
        }

        // dense output handling
        double hInt = ri_bs->maxStep;
        const GraggBulirschStoerStateInterpolator interpolator;
        if (! reject) {

            // extrapolate state at middle point of the step
            for (int j = 1; j <= k; ++j) {
                extrapolate(ri_bs, 0, j, diagonal, yMidDots[0], y_length);
            }

            const int mu = 2 * k - ri_bs->mudif + 3;

            for (int l = 0; l < mu; ++l) {

                // derivative at middle point of the step
                const int l2 = l / 2;
                double factor = pow(0.5 * ri_bs->sequence[l2], l);
                int fk_l2_length = ri_bs->sequence[l2] + 1;
                int middleIndex = fk_l2_length / 2;
                for (int i = 0; i < y_length; ++i) {
                    yMidDots[l + 1][i] = factor * fk[l2][middleIndex + l][i];
                }
                for (int j = 1; j <= k - l2; ++j) {
                    factor = pow(0.5 * ri_bs->sequence[j + l2], l);
                    int fk_l2j_length = ri_bs->sequence[l2+j] + 1;
                    middleIndex = fk_l2j_length / 2;
                    for (int i = 0; i < y_length; ++i) {
                        diagonal[j - 1][i] = factor * fk[l2 + j][middleIndex + l][i];
                    }
                    extrapolate(ri_bs, l2, j, diagonal, yMidDots[l + 1], y_length);
                }
                for (int i = 0; i < y_length; ++i) {
                    yMidDots[l + 1][i] *= ri_bs->stepSize;
                }

                // compute centered differences to evaluate next derivatives
                for (int j = (l + 1) / 2; j <= k; ++j) {
                    int fk_j_length = ri_bs->sequence[j] + 1;
                    for (int m = fk_j_length - 1; m >= 2 * (l + 1); --m) {
                        for (int i = 0; i < y_length; ++i) {
                            fk[j][m][i] -= fk[j][m - 2][i];
                        }
                    }
                }

            }

            // state at end of step
            const ODEStateAndDerivative stepEnd =
                equations.getMapper().mapStateAndDerivative(nextT, y1, computeDerivatives(nextT, y1));

            // set up interpolator covering the full step
            interpolator = new GraggBulirschStoerStateInterpolator(forward,
                    getStepStart(), stepEnd,
                    getStepStart(), stepEnd,
                    equations.getMapper(),
                    yMidDots, mu);

            if (mu >= 0 && ri_bs->useInterpolationError) {
                // use the interpolation error to limit stepsize
                const double interpError = interpolator.estimateError(scale);
                hInt = fabs(ri_bs->stepSize /
                        MAX(pow(interpError, 1.0 / (mu + 4)), 0.01));
                if (interpError > 10.0) {
                    hNew   = filterStep(ri_bs, hInt, forward, 0);
                    reject = 1;
                }
            }

        } else {
            interpolator = NULL;
        }

        if (! reject) {

            // Discrete events handling
            setStepStart(acceptStep(interpolator, finalTime));

            // prepare next step
            // beware that y1 is not always valid anymore here,
            // as some event may have triggered a reset
            // so we need to copy the new step start set previously
            y = getStepStart().getCompleteState();

            int optimalIter;
            if (k == 1) {
                optimalIter = 2;
                if (previousRejected) {
                    optimalIter = 1;
                }
            } else if (k <= targetIter) {
                optimalIter = k;
                if (ri_bs->costPerTimeUnit[k - 1] < ri_bs->orderControl1 * ri_bs->costPerTimeUnit[k]) {
                    optimalIter = k - 1;
                } else if (ri_bs->costPerTimeUnit[k] < ri_bs->orderControl2 * ri_bs->costPerTimeUnit[k - 1]) {
                    optimalIter = MIN(k + 1, ri_bs->sequence_length - 2);
                }
            } else {
                optimalIter = k - 1;
                if ((k > 2) && (ri_bs->costPerTimeUnit[k - 2] < ri_bs->orderControl1 * ri_bs->costPerTimeUnit[k - 1])) {
                    optimalIter = k - 2;
                }
                if (ri_bs->costPerTimeUnit[k] < ri_bs->orderControl2 * ri_bs->costPerTimeUnit[optimalIter]) {
                    optimalIter = MIN(k, ri_bs->sequence_length - 2);
                }
            }

            if (previousRejected) {
                // after a rejected step neither order nor stepsize
                // should increase
                targetIter = MIN(optimalIter, k);
                hNew = MIN(fabs(ri_bs->stepSize), ri_bs->optimalStep[targetIter]);
            } else {
                // stepsize control
                if (optimalIter <= k) {
                    hNew = filterStep(ri_bs, ri_bs->optimalStep[optimalIter], forward, 0);
                } else {
                    if ((k < targetIter) &&
                            (ri_bs->costPerTimeUnit[k] < ri_bs->orderControl2 * ri_bs->costPerTimeUnit[k - 1])) {
                        hNew = filterStep(ri_bs, ri_bs->optimalStep[k] * ri_bs->costPerStep[optimalIter + 1] / ri_bs->costPerStep[k], forward, 0);
                    } else {
                        hNew = filterStep(ri_bs, ri_bs->optimalStep[k] * ri_bs->costPerStep[optimalIter] / ri_bs->costPerStep[k], forward, 0);
                    }
                }

                targetIter = optimalIter;

            }

            newStep = 1;

        }

        hNew = MIN(hNew, hInt);
        if (! forward) {
            hNew = -hNew;
        }

        firstTime = 0;

        if (reject) {
            ri_bs->isLastStep = 0;
            previousRejected = 1;
        } else {
            previousRejected = 0;
        }

    } while (!isLastStep());

    const ODEStateAndDerivative finalState = getStepStart();
    resetInternalState();
    return finalState;

}



void reb_integrator_bs_part1(struct reb_simulation* r){
    r->t+=r->dt/2.;
}
void reb_integrator_bs_part2(struct reb_simulation* r){
    r->t+=r->dt/2.;
    r->dt_last_done = r->dt;
}

void reb_integrator_bs_synchronize(struct reb_simulation* r){
    // Do nothing.
}

void reb_integrator_bs_reset(struct reb_simulation* r){
    // Do nothing.
}
