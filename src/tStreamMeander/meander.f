c     "mndropt2.f"
C
c     version 5/12/97: convert program to subroutine called from 
c     'child'. 
c	changes:
c		stnserod--number of points in reach plus a 'tail'
c			of downstream points; replaces 'length'
c			as dimension of arrays
c		slope--now an array of slopes at each point
c		width-- "   "   "   "  channel widths "  "  "
c		nrough--  "
c		diam--    "
c		flow--    "
c		delta_x, delta_y, delta_z--displacement per time
c		grid variables--no longer needed
c               now calculated gaussian on the fly again.
c
c     version 2/4/97: various optimizations incl: calc. of 
c     'gaussian' array during initialization as 'look-up table'
c     for forcedist; use xs(*) vector in forcedist instead of
c     adding up dels's simplifies finding the downstream destination
c     point and the spreading.
c
c     version 3/06/96: modification of transverse slope formulation 
c     w.r.t. grain shear; "leftover" force discarded; channel erosion 
c     vector output instead of delx and dely; added spreaddeltas to 
c     "valleybounds," and latforce, lag, and spreaddeltas to
c     "rediscretize".
c
c     version 2/15/96: curvature calculated from delx and dely 
c     instead of phi; phi now unmodified, i.e., varies between -pi..pi.
c
c     stress induced by changing curvature is applied at some point 
c     downstream corresponding to the advection downstream during the 
c     time it takes for the induced stress to "move" crosstream.  the 
c     stress is also spread out according to a Gaussian, standard 
c     deviation constant; 
c
c     use Ikeda '89 for transverse slope calc.
c
c     version 1.6: eliminated erroneous division by dels in shear
c     stress calculation
c                 1.7  8/11: debugged version SL
c
c     $Id: meander.f,v 1.7 2003-05-09 17:06:07 childcvs Exp $
c
      subroutine meander (stations, stnserod, x, y, xs, dels, flow,  
     +                    rerody, lerody, slope, width, depth, 
     +                    diam, delta_x, delta_y, rightdepth, leftdepth,
     +                    lambda)
      implicit none
      integer stnserod, stations, dmnsn
      parameter (dmnsn = 6000)
      real*8 dels(*),
     +       phi(dmnsn), x(*), y(*), 
     +       delx(dmnsn), dely(dmnsn),  
     +       rho, grav, diam(*), flow(*), width(*), 
     +       slope(*), depth(*), vel(dmnsn), 
     +       curvature(dmnsn), 
     +       Acs(dmnsn), leftdepth(*), 
     +       rightdepth(*), deln(dmnsn), 
     +       latforce(dmnsn), lag(dmnsn), tauwall(dmnsn),
     +       spreaddelta_x(dmnsn), spreaddelta_y(dmnsn),
     +       xs(*), transfactor, transslope(dmnsn),
     +       rerody(*), lerody(*),
     +       delta_x(*), delta_y(*), lambda(*)
c-int intent(in) :: stations, stnserod, x, y, xs, dels, flow,  
c-int+     rerody, lerody, slope, diam
C
      continue
c     print *, 'stnserod in meander:', stnserod
c
C     declare parameter values; get initial channel config.:
c
      call initialize (stnserod, phi, x, y, delx, dely, 
     +                 rho, grav, width, lambda)
c     print *, 'stnserod in meander:', stnserod
c
C     define channel characteristics:
c
      call channel (stnserod, stations, flow, slope, diam, width, 
     +              rho, grav, phi, curvature, dels, 
     +              Acs, deln, rightdepth, leftdepth, depth, 
     +              vel, delx, dely, transfactor, transslope)
c
C     calculate lateral force and downstream lag:
c
      call forcelag (stnserod, stations, rho, vel, depth, width, 
     +               latforce, lag, curvature, Acs, dels, deln)
c
C     propagate and smooth lateral force:
c
      call forcedist (stnserod, stations, lambda, width, 
     +                lag, latforce, dels, phi, curvature, 
     +                depth, spreaddelta_x, spreaddelta_y,
     +                rightdepth, leftdepth, xs, tauwall)
c
C     change channel position:
c
      call changeposition (stnserod, stations, lerody, rerody,  
     +                     spreaddelta_x, spreaddelta_y, delx, dely, 
     +                     depth, delta_x, delta_y)
c      print *, 'last reach node K = ', transfactor
      return
      end
*
*
*
      subroutine initialize (stnserod, phi, x, y, 
     +                       delx, dely, rho, grav, width, lambda) 
      implicit none
      integer stnserod, s
      real*8 phi(*), x(*), y(*), 
     +       delx(*), dely(*), 
     +       rho, grav, width(*), lambda(*)
      real*8 angle
      external angle
c-int intent(in) :: stnserod, x, y
c
      continue
c     print *, 'stnserod in initialize:', stnserod
      rho = 1000.d0
      grav = 9.8d0
c      do s = 1, stnserod
c         lambda(s) = 2.75 * width(s)
c      end do
c     print *, 'stnserod in initialize:', stnserod
      do s = 1, stnserod - 1
c         print *, s, 'of', stnserod
         delx(s) = x(s + 1) - x(s)
         dely(s) = y(s + 1) - y(s)
         phi(s) = angle(dely(s), delx(s))
c         print *,s,delx(s),dely(s),phi(s),lambda(s)
      end do
      delx(stnserod) = delx(stnserod - 1)
      dely(stnserod) = dely(stnserod - 1)
      phi(stnserod) = phi(stnserod - 1)
c     print *, 'stnserod in initialize:', stnserod
      return
      end
*
*
*
      subroutine channel (stnserod, stations, flow, slope, diam, 
     +                    width, rho, grav, phi, 
     +                    curvature, dels, Acs, deln,
     +                    rightdepth, leftdepth, depth, vel, 
     +                    delx, dely, transfactor, transslope)
      implicit none
C
C  xmagtransslope is the absolute value of the transverse bedslope (transslope)
C  Acs is cross-sectional area of the half-channel at the inside of the bend
C
      integer stnserod, stations, s
      real*8 flow(*), slope(*), diam(*), width(*), rho, grav, 
     +       phi(*), curvature(*), 
     +       dels(*), Acs(*), deln(*),
     +       rightdepth(*), leftdepth(*), 
     +       transslope(*), depth(*), vel(*), shields, 
     +       critshields, 
     +       grainshields, corner, delx(*), 
     +       dely(*), transfactor, rectchan, radh, 
     +       xmagtransslope
      intrinsic abs, log10, log
c-int intent(in) :: stnserod, stations, flow, slope, diam, 
c-int+     width, rho, grav, phi, 
c-int+     dels,
c-int+     depth,
c-int+     delx, dely
      continue
C
      call getcurv (stnserod, stations, delx, dely, dels, curvature)
c     print *, 'stnserod in channel:', stnserod
      do s = 1, stnserod - 1
         if (slope(s) .le. 0.d0)
c     +        continue
     +       print *, 'neg. or zero slope:', slope(s), s, flow(s)
         if (width(s) .ne. 0.d0 .and. width(s) .ne. -2.d0*depth(s)) 
     +        then
c          radh = (width(s) * depth(s)) / (width(s) + 2.d0 * depth(s)) 
c          go back to H ~= R approx.
           radh = depth(s)
           if (depth(s) .le. 0.d0) then
              print *, depth(s)
              return
            end if
            vel(s) = flow(s) / depth(s) / width(s)
            shields = radh * slope(s) / 1.65d0 / diam(s) 
            if (diam(s) .lt. 0.00015d0) then
               critshields = 0.080d0 ! s*=1.01
            else if (diam(s) .lt. 0.00025d0) then
               critshields = 0.052d0 ! s*=2.84
            else if (diam(s) .lt. 0.00035d0) then
               critshields = 0.039d0 ! s*=5.22
            else if (diam(s) .lt. 0.00045d0) then
               critshields = 0.035d0 ! s*=8.04
            else if (diam(s) .lt. 0.00055d0) then
               critshields = 0.034d0 ! s*=11.2
            else if (diam(s) .lt. 0.00065d0) then
               critshields = 0.033d0 ! s*=14.8
            else if (diam(s) .lt. 0.00075d0) then
               critshields = 0.034d0 ! s*=18.6
            else if (diam(s) .lt. 0.00085d0) then
               critshields = 0.034d0 ! s*=22.7
            else if (diam(s) .lt. 0.00095d0) then
               critshields = 0.034d0 ! s*=27.1
            else if (diam(s) .lt. 0.0015d0) then
               critshields = 0.035d0 ! s*=31.8
            else if (diam(s) .lt. 0.0025d0) then
               critshields = 0.045d0 ! s*=89.9
            else if (diam(s) .lt. 0.0035d0) then
               critshields = 0.050d0 ! s*=165.
            else if (diam(s) .lt. 0.0045d0) then
               critshields = 0.055d0 ! s*=254.
            else
               critshields = 0.056d0 ! s*>=355
            end if
c           approximate Engelund diagram for subcritical flow:
            if (shields .ge. 0.1d0 .and. shields .lt. 1.d0 ) then
               grainshields = 10.d0 ** (0.74d0 
     +                        * (log10(shields) + 1.03d0) ** 2
     +                        - 1.18d0)
            else if (shields .ge. 1.d0 .and. shields .lt. 2.d0) then
               grainshields = 0.4d0 * shields ** 2
            else
               ! if (shields .ge. 2.0 .or. shields .lt. 0.1) then
               grainshields = shields
            end if
            corner = 2.d0 * depth(s) / width(s)
c          calculation for skin friction:
            transfactor = depth(s) * (grainshields / shields) 
     +                    * sqrt(grainshields / critshields) 
     +                    * (0.5695d0 * log(11.d0 * grainshields *  
     +                    depth(s) / shields / diam(s)) - 0.3606d0)
c           calculation for total shear:
c           transfactor = depth(s) * (shields / critshields) ** 0.5d0
c    +                    * (0.2279d0 / sqrt(grav * radh * slope(s)) * vel(s) 
c    +                    - 0.3606d0)
            rectchan = 0.5d0 * width(s) * depth(s)
            transslope(s) = transfactor * curvature(s) 
            xmagtransslope = abs(transslope(s))
            if (xmagtransslope .le. corner) then
               Acs(s) = rectchan - width(s) ** 2 * 
     +                  xmagtransslope / 8.d0
               rightdepth(s) = depth(s) + width(s) * 
     +                         transslope(s) / 2.d0
               leftdepth(s) = depth(s) - width(s) * transslope(s) / 2.d0
            else
               Acs(s) = 0.5d0 * depth(s) ** 2 / xmagtransslope
               if (curvature(s) .lt. 0.d0) then
                  rightdepth(s) = 0.d0
                  leftdepth(s) = 2.d0 * depth(s)
               else
                  rightdepth(s) = 2.d0 * depth(s)
                  leftdepth(s) = 0.d0
               end if
            end if
            if (curvature(s) .lt. 0.d0) then
               deln(s) = 0.5d0 * width(s) + Acs(s) * 2.d0 
     +                   / (depth(s) + rightdepth(s))
            else
               deln(s) = 0.5d0 * width(s) + Acs(s) * 2.d0 
     +                   / (depth(s) + leftdepth(s))
            end if
         end if
      end do
c      transfactor = transfactor / depth(s)
      return
      end
*
*
*
      subroutine getcurv (stnserod, stations, delx, dely, dels, 
     +                    curvature)
      implicit none
c     use law of cosines to find magnitude of angle, use cross 
c     product to find sign of curvature
c
      integer stnserod, stations, s
      real*8 delx(*), dely(*), dels(*), curvature(*), mag, sn, a, b, c,
     +       carg
c-int intent(in) :: stnserod, stations, delx, dely, dels
      intrinsic acos
      continue
c
      do s = 2, stnserod - 1
         a = dels(s) 
         b = dels(s - 1) 
         if( a .eq. 0.d0 .or. b .eq. 0.d0 ) then
            print *, 'dels(s) or dels(s-1) equals zero'
            stop
         end if
         c = sqrt((delx(s - 1) + delx(s)) * (delx(s - 1) + delx(s)) + 
     +            (dely(s - 1) + dely(s)) * (dely(s - 1) + dely(s)))
         carg = (a * a + b * b - c * c) / (2.d0 * a * b)
         if( carg .lt. -1.d0 ) then
c            if( carg .gt. -1.0001 ) carg = -1.d0
             carg = -1.d0
         end if
         if( carg .gt. 1.d0 ) then
c            if( carg .lt. 1.0001 ) carg = 1.d0
            carg = 1.d0
         end if
c            print *, 'arccosine argument out of bounds: ',carg,
c     +               ' at s= ', s,' of ',stnserod,'; a,b,c= ',a,b,c,
c     +               '; delx(s-1),delx(s),dely(s-1),dely(s): ',
c     +               delx(s-1),delx(s),dely(s-1),dely(s),
c     +               '; last curv: ',curvature(s-1)
c            stop
c         end if
         mag = (acos(-1.d0) - acos(carg)) / ((a + b) / 2.d0)
         sn = sign(1.d0, delx(s - 1) * dely(s) - dely(s - 1) * delx(s))
         curvature(s) = mag * sn
cDEBUGGING!
c         curvature(s) = sin(0.7*s)
c         if( curvature(s) .lt. 0 ) curvature(s)= -curvature(s)
c         print *,s,curvature(s),delx(s),dely(s)
            
c      This is tan(theta)!
c         curvature(s) = (dely(s) * delx(s - 1) - delx(s) * dely(s - 1)) 
c     +              / (delx(s) * delx(s - 1) + dely(s) * dely(s - 1))
c     +              / dels(s)
      end do
      curvature(1) = 0.d0
      curvature(stnserod) = 0.d0
c      curvature(stnserod - 1) = 0.d0
      return
      end
*
*
*
      subroutine forcelag (stnserod, stations, rho, vel, depth, width, 
     +                     latforce, lag, curvature, Acs, dels, deln)
      implicit none
c
C     calculate lateral force and downstream lag:
c
      integer stnserod, stations, s
      real*8 rho, vel(*), depth(*), width(*), latforce(*), 
     +       lag(*), curvature(*), 
     +       Acs(*), dels(*), deln(*), 
     +       forcefactor, xmagcurvep1, xmagcurve, latvel,
     +       delAcs
      intrinsic abs, cos
c-int intent(in) :: stnserod, stations, rho, vel, depth, width, 
c-int+     curvature, Acs, dels, deln
      continue
c
      if (stations .ne. stnserod) then
         do s = 1, stations - 1
            if (width(s) .ne. 0.d0 .and. 
     +          width(s) .gt. 2.d0 * depth(s)) then
               forcefactor = rho * vel(s) ** 2 / depth(s) 
     +              * (1.d0 - 2.d0 * depth(s) / width(s))
               latforce(s) = 0.d0
               lag(s) = 0.d0
               xmagcurvep1 = abs(curvature(s + 1)) 
               xmagcurve = abs(curvature(s)) 
               if ((curvature(s + 1) * curvature(s) .gt. 0.d0 .and.
     +           xmagcurvep1 .gt. xmagcurve .and.
     +           xmagcurvep1 * width(s) .le. 2.d0) .or.
     +          (curvature(s) .eq. 0.d0 .and.
     +           xmagcurvep1 * width(s) .le. 2.d0)) then
                  delAcs = (Acs(s + 1) - Acs(s))
               else if (curvature(s + 1) * curvature(s) .lt. 0.d0 .and.
     +                  xmagcurvep1 * width(s) .le. 2.d0) then
                  delAcs = Acs(s + 1) - depth(s) * width(s) / 2.d0
               else
                  delAcs = 0.d0
               end if
               latforce(s) = forcefactor * delAcs * delAcs 
     +                    / dels(s) * sign(1.d0, -curvature(s + 1)) 
     +                    * cos(dels(s) * xmagcurve / 2.d0)
               latvel = -1.d0 * vel(s) * delAcs / depth(s) / dels(s)
               if (latvel .gt. 0.d0) then
                     lag(s) = vel(s) * (deln(s + 1) + deln(s)) / 
     +                        2.d0 / latvel 
               else
                  latforce(s) = 0.d0
                  lag(s) = 0.d0
               end if
            else
               latforce(s) = 0.d0
               lag(s) = 0.d0
            end if
         end do
      else
         do s = 1, stations - 1
            if (width(s) .ne. 0.d0) then
               forcefactor = rho * vel(s) ** 2 / depth(s) 
     +              * (1.d0 - 2.d0 * depth(s) / width(s))
               latforce(s) = 0.d0
               lag(s) = 0.d0
               xmagcurvep1 = abs(curvature(s + 1)) 
               xmagcurve = abs(curvature(s)) 
               if ((curvature(s + 1) * curvature(s) .gt. 0.d0 .and.
     +           xmagcurvep1 .gt. xmagcurve .and.
     +           xmagcurvep1 * width(s) .le. 2.d0) .or.
     +          (curvature(s) .eq. 0.d0 .and.
     +           xmagcurvep1 * width(s) .le. 2.d0)) then
                  delAcs = (Acs(s + 1) - Acs(s))
               else if (curvature(s + 1) * curvature(s) .lt. 0.d0 .and.
     +               xmagcurvep1 * width(s) .le. 2.d0) then
                  delAcs = Acs(s + 1) - depth(s) * width(s) / 2.d0
               else
                  delAcs = 0.d0
               end if
               latforce(s) = forcefactor * delAcs * delAcs 
     +                       / dels(s) * sign(1.d0, -curvature(s + 1)) 
     +                    * cos(dels(s) * xmagcurve / 2.d0)
               latvel = -1.d0 * vel(s) * delAcs / depth(s) / dels(s)
               if (latvel .gt. 0.d0) then
                     lag(s) = vel(s) * (deln(s + 1) + deln(s)) / 
     +                        2.d0 / latvel 
               else
                  latforce(s) = 0.d0
                  lag(s) = 0.d0
               end if
            else
               latforce(s) = 0.d0
               lag(s) = 0.d0
             end if
         end do
      end if
      return
      end
*
*
*
      subroutine forcedist (stnserod, stations, lambda, width, 
     +                      lag, latforce, dels, phi, curvature, 
     +                      depth, spreaddelta_x, spreaddelta_y,
     +                      rightdepth, leftdepth, xs, tauwall)
      implicit none
C
      integer stnserod, stations, s, sp
      real*8 lambda(*), width(*), lag(*), 
     +       latforce(*), dels(*), phi(*), 
     +       curvature(*), depth(*), spreaddelta_x(*), 
     +       spreaddelta_y(*), rightdepth(*),
     +       leftdepth(*), xs(*),
     +       tauwall(*), gaussfactor, xstrt,
     +       xdel, xdepth, gaussian, xdest, tenlambda, xtrmnt
      intrinsic abs, cos, sin, exp
c-int intent(in) :: stnserod, stations, lambda, width, 
c-int+     lag, latforce, dels, phi, curvature, 
c-int+     depth,
c-int+     rightdepth, leftdepth, xs
      continue
c
      do s = 1, stnserod
         tauwall(s) = 0.d0
         spreaddelta_x(s) = 0.d0
         spreaddelta_y(s) = 0.d0
      end do
c     For each latforce(s), scan sp = s++ until the
c     distance downstream is greater than xs(s) + lag(s) 
c     + 2*lambda; for each sp, add to tauwall(sp): 
c        latforce(s) * gaussian(xs(s) + lag(s) - xs(sp))
c     find deeper bank; calc. 
c     normalized (unit vector divided by downstream
c     increment and bank depth) lateral direction vectors:
      do s = 1, stations
         tenlambda = 10.d0 * lambda(s)
         if (tenlambda .gt. xs(stations)) tenlambda = xs(stations)
         if (lag(s) .le. tenlambda) then
            xdest = xs(s) + lag(s)
            xtrmnt = xdest + 2.d0 * lambda(s)
            xstrt = xdest - 2.d0 * lambda(s)
            if (xstrt .lt. xs(s)) xstrt = xs(s)
            sp = s
            do while ( sp .le. stnserod )
               if (.not.( xs(sp) .le. xtrmnt)) exit
               if (xs(sp) .ge. xstrt) then
                  xdel = abs(xdest - xs(sp))
                  if (lambda(s) .ne. 0.d0) then
                     gaussian = exp(-1.d0 * xdel ** 2 / 2.d0 
     +                    / lambda(s) ** 2) / sqrt(2.d0 * 3.1416d0) 
     +                    / lambda(s) 
                     gaussfactor = gaussian * latforce(s)
                  else
                     gaussfactor = 0.d0
                  end if
                  tauwall(sp) = tauwall(sp) + gaussfactor
               end if
               sp = sp + 1
            end do
         end if
      end do
      do s = 1, stnserod
         xdepth = rightdepth(s)
         if (curvature(s) .lt. 0.d0) xdepth = leftdepth(s)
         if (xdepth .ne. 0.d0) then
            tauwall(s) = tauwall(s) / xdepth
            spreaddelta_x(s) = tauwall(s) * (-1.d0) * sin(phi(s))
            spreaddelta_y(s) = tauwall(s) * cos(phi(s))
         else
            tauwall(s) = 0.d0
            spreaddelta_x(s) = 0.d0
            spreaddelta_y(s) = 0.d0
         end if
c         print *,s,tauwall(s),spreaddelta_x(s),spreaddelta_y(s)
      end do
      return
      end
*
*
*
C
C  CHANGEPOSITION
C
C  This subroutine moves the meander nodes.
C  xcp is the cross-product of the flow direction and the direction of
C  bank erosion (perpendicular to the flow direction); the sign tells
C  you whether to use right or left bank erodibility.
C
      subroutine changeposition (stnserod, stations, lerody, rerody,  
     +                     spreaddelta_x, spreaddelta_y, delx, dely, 
     +                     depth, delta_x, delta_y)
      implicit none
c
C     change channel position:
c
      integer stnserod, stations, s
      real*8 lerody(*), rerody(*), spreaddelta_x(*), 
     +       spreaddelta_y(*), delx(*), dely(*), 
     +       depth(*), delta_x(*), delta_y(*), xcp
c-int intent(in) :: stnserod, stations, lerody, rerody,  
c-int+     spreaddelta_x, spreaddelta_y, delx, dely, depth
c
      continue
c      print *,'CHANGE CHANNEL POSITION:'
      do s = 1, stnserod
c        print *, 'ler ', lerody(s),'  rer',rerody(s)
         xcp = delx(s) * spreaddelta_y(s) - spreaddelta_x(s) * dely(s)
         if (xcp .lt. 0.d0) then
C            delta_x(s) = rerody(s) * depth(s) * spreaddelta_x(s) BUG FIX
C            delta_y(s) = rerody(s) * depth(s) * spreaddelta_y(s) 3/99 GT
            delta_x(s) = rerody(s) * spreaddelta_x(s)
            delta_y(s) = rerody(s) * spreaddelta_y(s)
         else
            delta_x(s) = lerody(s) * spreaddelta_x(s)
            delta_y(s) = lerody(s) * spreaddelta_y(s)
         end if
c         print *,s,delta_x(s),delta_y(s)
      end do
      return
      end
*
*
*
      real*8 function angle (y, x)
      implicit none
      real*8 y, x
      intrinsic atan2
c-int intent(in) :: y, x
      continue
      if (x .ne. 0.0d0 .or. y .ne. 0.0d0) then
         angle = atan2(y, x)
      else
         angle = 0.0d0
      end if
c      if (x .lt. 0.0) then
c         if (y .gt. 0.0) then
c            angle = 3.1416 + atan(y/x)
c         else
c            angle = -1.0 * (3.1416 - atan(y/x))
c         end if
c      else if (x .eq. 0.0) then
c         if (y .gt. 0.0) then
c            angle = 3.1416 / 2.0
c         else if (y .lt. 0.0) then
c            angle = -1.0 * 3.1416 / 2.0
c         else
c            angle = 0.0
c         end if
c      else
c         angle = atan(y/x)
c      end if
      return
      end
*
*
*


