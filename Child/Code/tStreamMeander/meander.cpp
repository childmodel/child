/**************************************************************************/
/**
**  @file
**  @brief definition of meander_
**
** translated by f2c (version 20020621).
**   Arnaud Desitter -- massaged somewhat to removed f2c dependencies.
*/
/**************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>

#include "meander.h"
#include "../Definitions.h"

#define integer int
#define doublereal double

static
void forcelag_(const integer *stnserod, const integer *stations, 
	       const doublereal *rho, const doublereal *vel,
	       const doublereal *depth, const doublereal *width,
	       doublereal *latforce, doublereal *lag,
	       const doublereal *curvature, const doublereal *acs,
	       const doublereal *dels, const doublereal *deln);
static
void forcedist_(const integer *stnserod, const integer *stations, 
		const doublereal *lambda, const doublereal *width,
		const doublereal *lag, const doublereal *latforce,
		const doublereal *dels, const doublereal *phi,
		const doublereal *curvature, const doublereal *depth,
		doublereal *spreaddelta_x__, doublereal *spreaddelta_y__,
		const doublereal *rightdepth, const doublereal *leftdepth, 
		const doublereal *xs, doublereal *tauwall);
static
void initialize_(const integer *stnserod, doublereal *phi, 
		 const doublereal *x, const doublereal *y,
		 doublereal *delx, doublereal *dely, 
		 doublereal *rho, doublereal *grav,
		 const doublereal *width, doublereal *lambda);
static
void changeposition_(const integer *stnserod, const integer * /*stations*/, 
		     const doublereal *lerody, const doublereal *rerody,
		     const doublereal *spreaddelta_x__,
		     const doublereal *spreaddelta_y__,
		     const doublereal *delx, const doublereal *dely,
		     const doublereal *depth,
		     doublereal *delta_x__, doublereal *delta_y__);
static
void channel_(const integer *stnserod, const integer *stations, 
	      const doublereal *flow, const doublereal *slope,
	      const doublereal *diam, const doublereal *width,
	      const doublereal * /*rho*/, const doublereal * /*grav*/,
	      const doublereal *phi, doublereal *curvature,
	      const doublereal *dels, 
	      doublereal *acs, doublereal *deln, doublereal *rightdepth,
	      doublereal *leftdepth, const doublereal *depth, 
	      doublereal *vel, const doublereal *delx, const doublereal *dely,
	      doublereal *transfactor, doublereal *transslope);
static
doublereal angle_(const doublereal *, const doublereal *);
static
void getcurv_(const integer *, const integer *, const doublereal *, 
	      const doublereal *, const doublereal *, doublereal *);


static
double d_sign(const doublereal *a, const doublereal *b){
  return *b<0 ? -fabs(*a) : fabs(*a);
}

/* Table of constant values */

static const doublereal c_b7 = 1.;

/*     "mndropt2.f" */

/*     version 6/9/2003 (SL): Updated code to 1999 version used in */
/*    Lancaster & Bras, 2002: */
/*       -curvature calculation had already been fixed; */
/*       -removed vertical flow momentum from force calc. (Jim Pizzuto */
/*        pointed out we shouldn't add perpendicular vectors!); */
/*       -removed multiplication of force by cosine of half the change */
/*        in angle (difficult to explain and physically justify); */
/*       -changed bank shear stress's depth divisor to the average */
/*        depth rather than the deeper bank depth (again, difficult */
/*        to explain and justify why we should always divide by the */
/*        "outer" bank depth; instead, just keep it simple and use */
/*        the avg. value); */

/*     version 5/12/97: convert program to subroutine called from */
/*     'child'. */
/* 	changes: */
/* 		stnserod--number of points in reach plus a 'tail' */
/* 			of downstream points; replaces 'length' */
/* 			as dimension of arrays */
/* 		slope--now an array of slopes at each point */
/* 		width-- "   "   "   "  channel widths "  "  " */
/* 		nrough--  " */
/* 		diam--    " */
/* 		flow--    " */
/* 		delta_x, delta_y, delta_z--displacement per time */
/* 		grid variables--no longer needed */
/*               now calculated gaussian on the fly again. */

/*     version 2/4/97: various optimizations incl: calc. of */
/*     'gaussian' array during initialization as 'look-up table' */
/*     for forcedist; use xs(*) vector in forcedist instead of */
/*     adding up dels's simplifies finding the downstream destination */
/*     point and the spreading. */

/*     version 3/06/96: modification of transverse slope formulation */
/*     w.r.t. grain shear; "leftover" force discarded; channel erosion */
/*     vector output instead of delx and dely; added spreaddeltas to */
/*     "valleybounds," and latforce, lag, and spreaddeltas to */
/*     "rediscretize". */

/*     version 2/15/96: curvature calculated from delx and dely */
/*     instead of phi; phi now unmodified, i.e., varies between -pi..pi. */

/*     stress induced by changing curvature is applied at some point */
/*     downstream corresponding to the advection downstream during the */
/*     time it takes for the induced stress to "move" crosstream.  the */
/*     stress is also spread out according to a Gaussian, standard */
/*     deviation constant; */

/*     use Ikeda '89 for transverse slope calc. */

/*     version 1.6: eliminated erroneous division by dels in shear */
/*     stress calculation */
/*                 1.7  8/11: debugged version SL */

/*     $Id: meander.cpp,v 1.17 2004-06-16 13:37:44 childcvs Exp $ */

void meander_(const integer *stations, const integer *stnserod, 
	      const doublereal *x, const doublereal *y,
	      const doublereal *xs, const doublereal *dels, 
	      const doublereal *flow, const doublereal *rerody,
	      const doublereal *lerody, const doublereal *slope,
	      const doublereal *width, const doublereal *depth, const doublereal *diam, 
	      doublereal *delta_x__, doublereal *delta_y__,
	      doublereal *rightdepth, doublereal *leftdepth,
	      doublereal *lambda)
{
    doublereal rho, transfactor, grav;

    doublereal
      *delx = new doublereal[*stnserod],
      *dely = new doublereal[*stnserod],
      *phi  = new doublereal[*stnserod],
      *tauwall = new doublereal[*stnserod],
      *spreaddelta_x__ = new doublereal[*stnserod],
      *spreaddelta_y__ = new doublereal[*stnserod],
      *curvature = new doublereal[*stnserod],
      *vel = new doublereal[*stnserod-1],
      *transslope = new doublereal[*stnserod-1],
      *acs = new doublereal[*stnserod],
      *deln = new doublereal[*stnserod],
      *latforce = new doublereal[*stations],
      *lag = new doublereal[*stations];


    memset(acs,0,*stnserod*sizeof(doublereal));
    memset(deln,0,*stnserod*sizeof(doublereal));
    memset(latforce,0,*stations*sizeof(doublereal));
    memset(lag,0,*stations*sizeof(doublereal));

    /* Parameter adjustments */
    --lambda;
    --leftdepth;
    --rightdepth;
    --delta_y__;
    --delta_x__;
    --diam;
    --depth;
    --width;
    --slope;
    --lerody;
    --rerody;
    --flow;
    --dels;
    --xs;
    --y;
    --x;

    /* Function Body */
/*     print *, 'stnserod in meander:', stnserod */

/*     declare parameter values; get initial channel config.: */

    initialize_(stnserod, phi, &x[1], &y[1], delx, dely, &rho, &grav, &width[
	    1], &lambda[1]);
/*     print *, 'stnserod in meander:', stnserod */

/*     define channel characteristics: */

    channel_(stnserod, stations, &flow[1], &slope[1], &diam[1], &width[1], &
	    rho, &grav, phi, curvature, &dels[1], acs, deln, &rightdepth[1], &
	    leftdepth[1], &depth[1], vel, delx, dely, &transfactor, 
	    transslope);

/*     calculate lateral force and downstream lag: */

    forcelag_(stnserod, stations, &rho, vel, &depth[1], &width[1], latforce, 
	    lag, curvature, acs, &dels[1], deln);

/*     propagate and smooth lateral force: */

    forcedist_(stnserod, stations, &lambda[1], &width[1], lag, latforce, &
	    dels[1], phi, curvature, &depth[1], spreaddelta_x__, 
	    spreaddelta_y__, &rightdepth[1], &leftdepth[1], &xs[1], tauwall);

/*     change channel position: */

    changeposition_(stnserod, stations, &lerody[1], &rerody[1], 
	    spreaddelta_x__, spreaddelta_y__, delx, dely, &depth[1], &
	    delta_x__[1], &delta_y__[1]);
/*      print *, 'last reach node K = ', transfactor */

    delete [] delx;
    delete [] dely;
    delete [] phi ;
    delete [] tauwall;
    delete [] spreaddelta_x__;
    delete [] spreaddelta_y__;
    delete [] curvature;
    delete [] vel;
    delete [] transslope;
    delete [] acs;
    delete [] deln;
    delete [] latforce;
    delete [] lag;

    return;
} /* meander_ */




void initialize_(const integer *stnserod, doublereal *phi, 
		 const doublereal *x, const doublereal *y,
		 doublereal *delx, doublereal *dely, 
		 doublereal *rho, doublereal *grav,
		 const doublereal * /*width*/, doublereal * /*lambda*/)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer s;

    /* Parameter adjustments */
    //--lambda;
    //--width;
    --dely;
    --delx;
    --y;
    --x;
    --phi;

    /* Function Body */
/*     print *, 'stnserod in initialize:', stnserod */
    *rho = RHO;
    *grav = GRAV;
/*      do s = 1, stnserod */
/*         lambda(s) = 2.75 * width(s) */
/*      end do */
/*     print *, 'stnserod in initialize:', stnserod */
    i__1 = *stnserod - 1;
    for (s = 1; s <= i__1; ++s) {
/*         print *, s, 'of', stnserod */
	delx[s] = x[s + 1] - x[s];
	dely[s] = y[s + 1] - y[s];
	phi[s] = angle_(&dely[s], &delx[s]);
/*         print *,s,delx(s),dely(s),phi(s),lambda(s) */
    }
    delx[*stnserod] = delx[*stnserod - 1];
    dely[*stnserod] = dely[*stnserod - 1];
    phi[*stnserod] = phi[*stnserod - 1];
/*     print *, 'stnserod in initialize:', stnserod */
    return;
} /* initialize_ */




void channel_(const integer *stnserod, const integer *stations, 
	      const doublereal *flow, const doublereal *slope,
	      const doublereal *diam, const doublereal *width,
	      const doublereal * /*rho*/, const doublereal * /*grav*/,
	      const doublereal * /*phi*/, doublereal *curvature,
	      const doublereal *dels, 
	      doublereal *acs, doublereal *deln, doublereal *rightdepth,
	      doublereal *leftdepth, const doublereal *depth, 
	      doublereal *vel, const doublereal *delx, const doublereal *dely,
	      doublereal *transfactor, doublereal *transslope)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    doublereal rectchan;
    integer s;
    doublereal critshields, radh, grainshields, corner, xmagtransslope, 
	    shields;

/*  xmagtransslope is the absolute value of the transverse bedslope (transslope) */
/*  Acs is cross-sectional area of the half-channel at the inside of the bend */

    /* Parameter adjustments */
    --transslope;
    --dely;
    --delx;
    --vel;
    --depth;
    --leftdepth;
    --rightdepth;
    --deln;
    --acs;
    --dels;
    --curvature;
    //--phi;
    --width;
    --diam;
    --slope;
    --flow;

    /* Function Body */

    getcurv_(stnserod, stations, &delx[1], &dely[1], &dels[1], &curvature[1]);
/*     print *, 'stnserod in channel:', stnserod */
    i__1 = *stnserod - 1;
    for (s = 1; s <= i__1; ++s) {
	if (slope[s] <= 0.) {
	  std::cout << "neg. or zero slope:" << slope[s] << " " << s 
	       << " " << flow[s] << std::endl;
	}
	if (width[s] != 0. && width[s] != depth[s] * -2.) {
/*          radh = (width(s) * depth(s)) / (width(s) + 2.d0 * depth(s)) */
/*          go back to H ~= R approx. */
	    radh = depth[s];
	    if (depth[s] <= 0.) {
	      std::cout << depth[s] << std::endl;
		return;
	    }
	    vel[s] = flow[s] / depth[s] / width[s];
	    shields = radh * slope[s] / 1.65 / diam[s];
	    if (diam[s] < 1.5e-4) {
		critshields = .08;
/* s*=1.01 */
	    } else if (diam[s] < 2.5e-4) {
		critshields = .052;
/* s*=2.84 */
	    } else if (diam[s] < 3.5e-4) {
		critshields = .039;
/* s*=5.22 */
	    } else if (diam[s] < 4.5e-4) {
		critshields = .035;
/* s*=8.04 */
	    } else if (diam[s] < 5.5e-4) {
		critshields = .034;
/* s*=11.2 */
	    } else if (diam[s] < 6.5e-4) {
		critshields = .033;
/* s*=14.8 */
	    } else if (diam[s] < 7.5e-4) {
		critshields = .034;
/* s*=18.6 */
	    } else if (diam[s] < 8.5e-4) {
		critshields = .034;
/* s*=22.7 */
	    } else if (diam[s] < 9.5e-4) {
		critshields = .034;
/* s*=27.1 */
	    } else if (diam[s] < .0015) {
		critshields = .035;
/* s*=31.8 */
	    } else if (diam[s] < .0025) {
		critshields = .045;
/* s*=89.9 */
	    } else if (diam[s] < .0035) {
		critshields = .05;
/* s*=165. */
	    } else if (diam[s] < .0045) {
		critshields = .055;
/* s*=254. */
	    } else {
		critshields = .056;
/* s*>=355 */
	    }
/*           approximate Engelund diagram for subcritical flow: */
	    if (shields >= .1 && shields < 1.) {
/* Computing 2nd power */
		d__2 = log10(shields) + 1.03;
		d__1 = d__2 * d__2 * .74 - 1.18;
		grainshields = pow(10., d__1);
	    } else if (shields >= 1. && shields < 2.) {
/* Computing 2nd power */
		d__1 = shields;
		grainshields = d__1 * d__1 * .4;
	    } else {
/* if (shields .ge. 2.0 .or. shields .lt. 0.1) then */
		grainshields = shields;
	    }
	    corner = depth[s] * 2. / width[s];
/*          calculation for skin friction: */
	    *transfactor = depth[s] * (grainshields / shields) * sqrt(
		    grainshields / critshields) * (log(grainshields * 11. * 
		    depth[s] / shields / diam[s]) * .5695 - .3606);
/*           calculation for total shear: */
/*           transfactor = depth(s) * (shields / critshields) ** 0.5d0 */
/*    +                    * (0.2279d0 / sqrt(grav * radh * slope(s)) * vel(s) */
/*    +                    - 0.3606d0) */
	    rectchan = width[s] * .5 * depth[s];
	    transslope[s] = *transfactor * curvature[s];
	    xmagtransslope = fabs(transslope[s]);
	    if (xmagtransslope <= corner) {
/* Computing 2nd power */
		d__1 = width[s];
		acs[s] = rectchan - d__1 * d__1 * xmagtransslope / 8.;
		rightdepth[s] = depth[s] + width[s] * transslope[s] / 2.;
		leftdepth[s] = depth[s] - width[s] * transslope[s] / 2.;
	    } else {
/* Computing 2nd power */
		d__1 = depth[s];
		acs[s] = d__1 * d__1 * .5 / xmagtransslope;
		if (curvature[s] < 0.) {
		    rightdepth[s] = 0.;
		    leftdepth[s] = depth[s] * 2.;
		} else {
		    rightdepth[s] = depth[s] * 2.;
		    leftdepth[s] = 0.;
		}
	    }
	    if (curvature[s] < 0.) {
		deln[s] = width[s] * .5 + acs[s] * 2. / (depth[s] + 
			rightdepth[s]);
	    } else {
		deln[s] = width[s] * .5 + acs[s] * 2. / (depth[s] + leftdepth[
			s]);
	    }
	}
    }
/*      transfactor = transfactor / depth(s) */
    return;
} /* channel_ */




void getcurv_(const integer *stnserod, const integer * /*stations*/, 
	      const doublereal *delx, const doublereal *dely,
	      const doublereal *dels, doublereal *curvature)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    doublereal a, b, c__;
    integer s;
    doublereal sn, mag, carg;

/*     use law of cosines to find magnitude of angle, use cross */
/*     product to find sign of curvature */

    /* Parameter adjustments */
    --curvature;
    --dels;
    --dely;
    --delx;

    /* Function Body */

    i__1 = *stnserod - 1;
    for (s = 2; s <= i__1; ++s) {
	a = dels[s];
	b = dels[s - 1];
	if (a == 0. || b == 0.) {
	  std::cout << "dels(s) or dels(s-1) equals zero" << std::endl;
	  exit(1);
	}
	c__ = sqrt((delx[s - 1] + delx[s]) * (delx[s - 1] + delx[s]) + (dely[
		s - 1] + dely[s]) * (dely[s - 1] + dely[s]));
	carg = (a * a + b * b - c__ * c__) / (a * 2. * b);
	if (carg < -1.) {
/*            if( carg .gt. -1.0001 ) carg = -1.d0 */
	    carg = -1.;
	}
	if (carg > 1.) {
/*            if( carg .lt. 1.0001 ) carg = 1.d0 */
	    carg = 1.;
	}
/*            print *, 'arccosine argument out of bounds: ',carg, */
/*     +               ' at s= ', s,' of ',stnserod,'; a,b,c= ',a,b,c, */
/*     +               '; delx(s-1),delx(s),dely(s-1),dely(s): ', */
/*     +               delx(s-1),delx(s),dely(s-1),dely(s), */
/*     +               '; last curv: ',curvature(s-1) */
/*            stop */
/*         end if */
	mag = (acos(-1.) - acos(carg)) / ((a + b) / 2.);
	d__1 = delx[s - 1] * dely[s] - dely[s - 1] * delx[s];
	sn = d_sign(&c_b7, &d__1);
	curvature[s] = mag * sn;
/* DEBUGGING! */
/*         curvature(s) = sin(0.7*s) */
/*         if( curvature(s) .lt. 0 ) curvature(s)= -curvature(s) */
/*         print *,s,curvature(s),delx(s),dely(s) */
/*      This is tan(theta)! */
/*         curvature(s) = (dely(s) * delx(s - 1) - delx(s) * dely(s - 1)) */
/*     +              / (delx(s) * delx(s - 1) + dely(s) * dely(s - 1)) */
/*     +              / dels(s) */
    }
    curvature[1] = 0.;
    curvature[*stnserod] = 0.;
/*      curvature(stnserod - 1) = 0.d0 */
    return;
} /* getcurv_ */




void forcelag_(const integer *stnserod, const integer *stations, 
	       const doublereal *rho, const doublereal *vel,
	       const doublereal *depth, const doublereal *width,
	       doublereal *latforce, doublereal *lag,
	       const doublereal *curvature, const doublereal *acs,
	       const doublereal *dels, const doublereal *deln)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    integer s;
    doublereal xmagcurve, xmagcurvep1, forcefactor, delacs, latvel;


/*     calculate lateral force and downstream lag: */

    /* Parameter adjustments */
    --deln;
    --dels;
    --acs;
    --curvature;
    --lag;
    --latforce;
    --width;
    --depth;
    --vel;

    /* Function Body */

    if (*stations != *stnserod) {
	i__1 = *stations - 1;
	for (s = 1; s <= i__1; ++s) {
	    if (width[s] != 0. && width[s] > depth[s] * 2.) {
/* Computing 2nd power */
		d__1 = vel[s];
		forcefactor = *rho * (d__1 * d__1) / depth[s];
/* * (1. - depth[s] * 2. / width[s]);*/
		latforce[s] = 0.;
		lag[s] = 0.;
		xmagcurvep1 = fabs(curvature[s + 1]);
		xmagcurve =fabs(curvature[s]);
		if ((curvature[s + 1] * curvature[s] > 0. && xmagcurvep1 > 
			xmagcurve && xmagcurvep1 * width[s] <= 2.) || 
			(curvature[s] == 0. && xmagcurvep1 * width[s] <= 2.)) {
		    delacs = acs[s + 1] - acs[s];
		} else if (curvature[s + 1] * curvature[s] < 0. && 
			xmagcurvep1 * width[s] <= 2.) {
		    delacs = acs[s + 1] - depth[s] * width[s] / 2.;
		} else {
		    delacs = 0.;
		}
		d__1 = -curvature[s + 1];
		latforce[s] = forcefactor * delacs * delacs / dels[s] * 
        d_sign(&c_b7, &d__1); /* * cos(dels[s] * xmagcurve / 2.); */
		latvel = vel[s] * -1. * delacs / depth[s] / dels[s];
		if (latvel > 0.) {
		    lag[s] = vel[s] * (deln[s + 1] + deln[s]) / 2. / latvel;
		} else {
		    latforce[s] = 0.;
		    lag[s] = 0.;
		}
	    } else {
		latforce[s] = 0.;
		lag[s] = 0.;
	    }
	}
    } else {
	i__1 = *stations - 1;
	for (s = 1; s <= i__1; ++s) {
	    if (width[s] != 0.) {
/* Computing 2nd power */
		d__1 = vel[s];
		forcefactor = *rho * (d__1 * d__1) / depth[s];
/* * (1. - depth[s] * 2. / width[s]); */
		latforce[s] = 0.;
		lag[s] = 0.;
		xmagcurvep1 = fabs(curvature[s + 1]);
		xmagcurve = fabs(curvature[s]);
		if ((curvature[s + 1] * curvature[s] > 0. && xmagcurvep1 > 
			xmagcurve && xmagcurvep1 * width[s] <= 2.) || 
			(curvature[s] == 0. && xmagcurvep1 * width[s] <= 2.)) {
		    delacs = acs[s + 1] - acs[s];
		} else if (curvature[s + 1] * curvature[s] < 0. && 
			xmagcurvep1 * width[s] <= 2.) {
		    delacs = acs[s + 1] - depth[s] * width[s] / 2.;
		} else {
		    delacs = 0.;
		}
		d__1 = -curvature[s + 1];
		latforce[s] = forcefactor * delacs * delacs / dels[s] * 
        d_sign(&c_b7, &d__1); /* * cos(dels[s] * xmagcurve / 2.); */
		latvel = vel[s] * -1. * delacs / depth[s] / dels[s];
		if (latvel > 0.) {
		    lag[s] = vel[s] * (deln[s + 1] + deln[s]) / 2. / latvel;
		} else {
		    latforce[s] = 0.;
		    lag[s] = 0.;
		}
	    } else {
		latforce[s] = 0.;
		lag[s] = 0.;
	    }
	}
    }
    return;
} /* forcelag_ */




void forcedist_(const integer *stnserod, const integer *stations, 
		const doublereal *lambda, const doublereal * /*width*/,
		const doublereal *lag, const doublereal *latforce,
		const doublereal * /*dels*/, const doublereal *phi,
		const doublereal *curvature, const doublereal *depth,
		doublereal *spreaddelta_x__, doublereal *spreaddelta_y__,
		const doublereal *rightdepth, const doublereal *leftdepth, 
		const doublereal *xs, doublereal *tauwall)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    doublereal gaussian;
    integer s;
    doublereal tenlambda;
    integer sp;
    doublereal gaussfactor, xdel, xdest, xstrt, /*xdepth,*/ xtrmnt;


    /* Parameter adjustments */
    --tauwall;
    --xs;
    --leftdepth;
    --rightdepth;
    --spreaddelta_y__;
    --spreaddelta_x__;
    --depth;
    --curvature;
    --phi;
    //--dels;
    --latforce;
    --lag;
    //--width;
    --lambda;

    /* Function Body */

    i__1 = *stnserod;
    for (s = 1; s <= i__1; ++s) {
	tauwall[s] = 0.;
	spreaddelta_x__[s] = 0.;
	spreaddelta_y__[s] = 0.;
    }
/*     For each latforce(s), scan sp = s++ until the */
/*     distance downstream is greater than xs(s) + lag(s) */
/*     + 2*lambda; for each sp, add to tauwall(sp): */
/*        latforce(s) * gaussian(xs(s) + lag(s) - xs(sp)) */
/*     find deeper bank; calc. */
/*     normalized (unit vector divided by downstream */
/*     increment and bank depth) lateral direction vectors: */
    i__1 = *stations;
    for (s = 1; s <= i__1; ++s) {
	tenlambda = lambda[s] * 10.;
	if (tenlambda > xs[*stations]) {
	    tenlambda = xs[*stations];
	}
	if (lag[s] <= tenlambda) {
	    xdest = xs[s] + lag[s];
	    xtrmnt = xdest + lambda[s] * 2.;
	    xstrt = xdest - lambda[s] * 2.;
	    if (xstrt < xs[s]) {
		xstrt = xs[s];
	    }
	    sp = s;
	    while( sp <= *stnserod && xs[sp] <= xtrmnt) {
		if (xs[sp] >= xstrt) {
		    xdel = fabs(xdest - xs[sp]);
		    if (lambda[s] != 0.) {
/* Computing 2nd power */
			d__1 = xdel;
/* Computing 2nd power */
			d__2 = lambda[s];
			gaussian = exp(d__1 * d__1 * -1. / 2. / (d__2 * d__2))
				 / sqrt(2. * 3.1416) / lambda[s];
			gaussfactor = gaussian * latforce[s];
		    } else {
			gaussfactor = 0.;
		    }
		    tauwall[sp] += gaussfactor;
		}
		++sp;
	    }
	}
    }
    i__1 = *stnserod;
    for (s = 1; s <= i__1; ++s) {
       /*xdepth = rightdepth[s];*/
       /*if (curvature[s] < 0.) {*/
       /*xdepth = leftdepth[s];*/
       /*}*/
       /*if (xdepth != 0.) */
       if(depth[s] != 0.){
          /*tauwall[s] /= xdepth;*/
          tauwall[s] /= depth[s];
          spreaddelta_x__[s] = tauwall[s] * -1. * sin(phi[s]);
          spreaddelta_y__[s] = tauwall[s] * cos(phi[s]);
       } else {
          tauwall[s] = 0.;
          spreaddelta_x__[s] = 0.;
          spreaddelta_y__[s] = 0.;
       }
/*         print *,s,tauwall(s),spreaddelta_x(s),spreaddelta_y(s) */
    }
    return;
} /* forcedist_ */





/*  CHANGEPOSITION */

/*  This subroutine moves the meander nodes. */
/*  xcp is the cross-product of the flow direction and the direction of */
/*  bank erosion (perpendicular to the flow direction); the sign tells */
/*  you whether to use right or left bank erodibility. */

void changeposition_(const integer *stnserod, const integer * /*stations*/, 
		     const doublereal *lerody, const doublereal *rerody,
		     const doublereal *spreaddelta_x__,
		     const doublereal *spreaddelta_y__,
		     const doublereal *delx, const doublereal *dely,
		     const doublereal * /*depth*/,
		     doublereal *delta_x__, doublereal *delta_y__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer s;
    doublereal xcp;


/*     change channel position: */


    /* Parameter adjustments */
    --delta_y__;
    --delta_x__;
    //--depth;
    --dely;
    --delx;
    --spreaddelta_y__;
    --spreaddelta_x__;
    --rerody;
    --lerody;

    /* Function Body */
/*      print *,'CHANGE CHANNEL POSITION:' */
    i__1 = *stnserod;
    for (s = 1; s <= i__1; ++s) {
/*        print *, 'ler ', lerody(s),'  rer',rerody(s) */
	xcp = delx[s] * spreaddelta_y__[s] - spreaddelta_x__[s] * dely[s];
	if (xcp < 0.) {
/*            delta_x(s) = rerody(s) * depth(s) * spreaddelta_x(s) BUG FIX */
/*            delta_y(s) = rerody(s) * depth(s) * spreaddelta_y(s) 3/99 GT */
	    delta_x__[s] = rerody[s] * spreaddelta_x__[s];
	    delta_y__[s] = rerody[s] * spreaddelta_y__[s];
	} else {
	    delta_x__[s] = lerody[s] * spreaddelta_x__[s];
	    delta_y__[s] = lerody[s] * spreaddelta_y__[s];
	}
/*         print *,s,delta_x(s),delta_y(s) */
    }
    return;
} /* changeposition_ */




doublereal angle_(const doublereal *y, const doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    if (*x != 0. || *y != 0.) {
	ret_val = atan2(*y, *x);
    } else {
	ret_val = 0.;
    }
/*      if (x .lt. 0.0) then */
/*         if (y .gt. 0.0) then */
/*            angle = 3.1416 + atan(y/x) */
/*         else */
/*            angle = -1.0 * (3.1416 - atan(y/x)) */
/*         end if */
/*      else if (x .eq. 0.0) then */
/*         if (y .gt. 0.0) then */
/*            angle = 3.1416 / 2.0 */
/*         else if (y .lt. 0.0) then */
/*            angle = -1.0 * 3.1416 / 2.0 */
/*         else */
/*            angle = 0.0 */
/*         end if */
/*      else */
/*         angle = atan(y/x) */
/*      end if */
    return ret_val;
} /* angle_ */

