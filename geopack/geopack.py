import numpy as np
from geopack import t89,t96,t01,t04
import os.path
import datetime

def init_igrf():
    """
    Initialize the IGRF coefficients and related coefs.
    Should be called once and only once when importing the geopack module.
    """

    global igrf, nmn,mns, nyear,years,yruts

    print('Load IGRF coefficients ...')

    bfn = 'igrf12coeffs.txt'
    locffn = os.path.join(os.path.dirname(__file__), bfn)

    nheader = 3
    with open(locffn, 'r') as file:
        for i in range(nheader):
            next(file)
        header = file.readline().rstrip()
        cols = header.split()
        # first 3 columns are g/h flag, n, m.
        years = np.array([np.int32(j[0:4]) for j in cols[3:]])
        nyear = len(years)
        lines = file.read().splitlines()

    cols = lines[-1].split()
    k = np.int32(cols[1]) + 1
    nmn = np.int32((k + 1) * k * 0.5)
    igrf = np.zeros((nmn, nyear, 2), dtype=float)
    mns = np.empty((nmn, 2), dtype=np.int32)
    l = 0
    for i in range(k):
        for j in range(i + 1):
            mns[l, :] = [i, j]
            l += 1

    for line in lines:
        cols = line.split()
        if cols[0] == 'g':
            i = 0
        else:
            i = 1
        n, m = np.int32(cols[1:3])
        mn = np.int32(n * (n + 1) * 0.5 + m)
        igrf[mn, :, i] = [np.float(j) for j in cols[3:]]

    # treat the last column
    years[-1] += 5
    igrf[:, -1, :] = igrf[:, -2, :] + igrf[:, -1, :] * 5
    yruts = np.empty(nyear, dtype=float)
    t0_datetime = datetime.datetime(1970,1,1)
    for i in range(nyear):
        yruts[i] = (datetime.datetime(years[i],1,1)-t0_datetime).total_seconds()
    # separate g/h
    igrf = {'g': igrf[:,:,0], 'h': igrf[:,:,1]}


def load_igrf(ut):
    """
    Load the IGRF coefficients for the given time.
    :param ut: ut sec, a float.
    :return: g,h. The IGRF coef at given time.
    """

    # igrf should be initilized already.
    global igrf, nmn,nyear, mns, years, yruts
    try: isinstance(igrf, dict)
    except: init_igrf()

    # locate the two years of interest.
    if ut <= yruts[0]: yridx = 0
    elif ut >= yruts[-1]: yridx = -2
    else:
        yridx = np.argwhere(yruts <= ut)[-1]

    g0 = np.squeeze(igrf['g'][:,yridx])
    h0 = np.squeeze(igrf['h'][:,yridx])
    g1 = np.squeeze(igrf['g'][:,yridx+1])
    h1 = np.squeeze(igrf['h'][:,yridx+1])

    ut0 = yruts[yridx]
    ut1 = yruts[yridx+1]
    f1 = (ut-ut0)/(ut1-ut0)
    f0 = 1-f1

    return g0*f0+g1*f1, h0*f0+h1*f1




def igrf_gsm(xgsm,ygsm,zgsm):
    """
    Calculates components of the main (internal) geomagnetic field in the geocentric solar
    magnetospheric coordinate system, using IAGA international geomagnetic reference model
    coefficients (e.g., http://www.ngdc.noaa.gov/iaga/vmod/igrf.html revised: 22 march, 2005)

    Before the first call of this subroutine, or if the date/time
    was changed, the model coefficients and GEO-GSM rotation matrix elements should be updated
    by calling the subroutine recalc

    Python version by Sheng Tian

    :param xgsm,ygsm,zgsm: cartesian GSM coordinates (in units Re=6371.2 km)
    :return: hxgsm,hygsm,hzgsm. Cartesian GSM components of the main geomagnetic field in nanotesla
    """

    xgeo,ygeo,zgeo = geogsm(xgsm,ygsm,zgsm, -1)
    r,theta,phi = sphcar(xgeo,ygeo,zgeo, -1)
    br,btheta,bphi = igrf_geo(r,theta,phi)
    bxgeo,bygeo,bzgeo = bspcar(theta,phi,br,btheta,bphi)
    return geogsm(bxgeo,bygeo,bzgeo, 1)


def igrf_geo(r,theta,phi):
    """
    Calculates components of the main (internal) geomagnetic field in the spherical geographic
    (geocentric) coordinate system, using IAGA international geomagnetic reference model
    coefficients  (e.g., http://www.ngdc.noaa.gov/iaga/vmod/igrf.html, revised: 22 march, 2005)

    Before the first call of this subroutine, or if the time was changed,
    the model coefficients should be updated by calling the subroutine recalc

    Python version by Sheng Tian

    :param r: spherical geographic (geocentric) coordinates: radial distance r in units Re=6371.2 km
    :param theta: colatitude theta in radians
    :param phi: longitude phi in radians
    :return: br, btheta, bphi. Spherical components of the main geomagnetic field in nanotesla
        (positive br outward, btheta southward, bphi eastward)
    """


    # common /geopack2/ g(105),h(105),rec(105)
    global g, h, rec

    ct = np.cos(theta)
    st = np.sin(theta)
    minst = 1e-5
    if np.abs(st) < minst: smlst = True
    else: smlst = False

    # In this new version, the optimal value of the parameter nm (maximal order of the spherical
    # harmonic expansion) is not user-prescribed, but calculated inside the subroutine, based
    # on the value of the radial distance r:
    irp3 = np.int64(r+2)
    nm = np.int64(3+30/irp3)
    if nm > 13: nm = 13
    k = nm+1

    # r dependence is encapsulated here.
    a = np.empty(k)
    b = np.empty(k)
    ar = 1/r        # a/r
    a[0] = ar*ar    # a[n] = (a/r)^(n+2).
    b[0] = a[0]     # b[n] = (n+1)(a/r)^(n+2)
    for n in range(1,k):
        a[n] = a[n-1]*ar
        b[n] = a[n]*(n+1)


    # t - short for theta, f - short for phi.
    br,bt,bf = [0.]*3
    d,p = [0.,1]

    # m = 0. P^n,0
    m = 0
    smf,cmf = [0.,1]
    p1,d1,p2,d2 = [p,d,0.,0]
    l0 = 0
    mn = l0
    for n in range(m,k):
        w = g[mn]*cmf+h[mn]*smf
        br += b[n]*w*p1          # p1 is P^n,m.
        bt -= a[n]*w*d1          # d1 is dP^n,m/dt.
        xk = rec[mn]
        # Eq 16c and its derivative on theta.
        d0 = ct*d1-st*p1-xk*d2   # dP^n,m/dt = ct*dP^n-1,m/dt - st*P_n-1,m - K^n,m*dP^n-2,m/dt
        p0 = ct*p1-xk*p2         # P^n,m = ct*P^n-1,m - K^n,m*P^n-2,m
        d2,p2,d1 = [d1,p1,d0]
        p1 = p0
        mn += n+1

    # Eq 16b and its derivative on theta.
    d = st*d+ct*p   # dP^m,m/dt = st*dP^m-1,m-1/dt + ct*P^m-1,m-1
    p = st*p        # P^m,m = st*P^m-1,m-1

    # Similarly for P^n,m
    l0 = 0
    for m in range(1,k):        # sum over m
        smf = np.sin(m*phi)     # sin(m*phi)
        cmf = np.cos(m*phi)     # cos(m*phi)
        p1,d1,p2,d2 = [p,d,0.,0]
        tbf = 0.
        l0 += m+1
        mn = l0
        for n in range(m,k):    # sum over n
            w=g[mn]*cmf+h[mn]*smf   # [g^n,m*cos(m*phi)+h^n,m*sin(m*phi)]
            br += b[n]*w*p1
            bt -= a[n]*w*d1
            tp = p1
            if smlst: tp = d1
            tbf += a[n]*(g[mn]*smf-h[mn]*cmf)*tp
            xk = rec[mn]
            d0 = ct*d1-st*p1-xk*d2   # dP^n,m/dt = ct*dP^n-1,m/dt - st*P_n-1,m - K^n,m*dP^n-2,m/dt
            p0 = ct*p1-xk*p2         # P^n,m = ct*P^n-1,m - K^n,m*P^n-2,m
            d2,p2,d1 = [d1,p1,d0]
            p1=p0
            mn += n+1

        d = st*d+ct*p
        p = st*p

        # update B_phi.
        tbf *= m
        bf += tbf

    if smlst:
        if ct < 0.: bf = -bf
    else: bf /= st

    return br,bt,bf


def dip(xgsm,ygsm,zgsm):
    """
    Calculates gsm components of a geodipole field with the dipole moment
    corresponding to the epoch, specified by calling subroutine recalc (should be
    invoked before the first use of this one and in case the date/time was changed).

    :param xgsm,ygsm,zgsm: GSM coordinates in Re (1 Re = 6371.2 km)
    :return: bxgsm,bygsm,gzgsm. Field components in gsm system, in nanotesla.

    Last modification: May 4, 2005.
    Author: N. A. Tsyganenko
    """

    # common /geopack1/ aaa(10),sps,cps,bbb(23)
    # common /geopack2/ g(105),h(105),rec(105)
    global aaa, sps,cps, bbb, g, h, rec

    dipmom = np.sqrt(g[1]**2+g[2]**2+h[2]**2)

    p = xgsm**2
    u = zgsm**2
    v = 3*zgsm*xgsm
    t = ygsm**2
    q = dipmom/np.sqrt(p+t+u)**5

    bxgsm = q*((t+u-2.*p)*sps-v*cps)
    bygsm = -3.*ygsm*q*(xgsm*sps+zgsm*cps)
    bzgsm = q*((p+t-2.*u)*cps-v*sps)

    return bxgsm,bygsm,bzgsm



def recalc(ut):
    """
    1. Prepares elements of rotation matrices for transformations of vectors between
        several coordinate systems, most frequently used in space physics.
    2. Prepares coefficients used in the calculation of the main geomagnetic field (igrf model)

    This subroutine should be invoked before using the following subroutines:
        igrf_geo, igrf_gsm, dip, geomag, geogsm, magsm, smgsm, gsmgse, geigeo.
    There is no need to repeatedly invoke recalc, if multiple calculations are made for the same date and time.

    :param ut: Universal time in second.
    :return: psi. Dipole tilt angle in radian.

    Python version by Sheng Tian
    """


    # The common block /geopack1/ contains elements of the rotation matrices and other
    # parameters related to the coordinate transformations performed by this package
    # common /geopack1/ st0,ct0,sl0,cl0,ctcl,stcl,ctsl,stsl,sfi,cfi,sps,
    # cps,shi,chi,hi,psi,xmut,a11,a21,a31,a12,a22,a32,a13,a23,a33,ds3,
    # cgst,sgst,ba(6)

    # st0/ct0 - sin/cos of teta0 (colat in geo?). sl0/cl0 - sin/cos of lambda0 (longitude in geo?).
    # ctcl/stcl - . ctsl/stsl - . geo x and z?
    # sfi/cfi/xmut - rotate angle between mag and sm and its sin/cos.
    # sps/cps/psi - tilt angle and its sin/cos.
    # shi/chi/hi - rotate angle between gse to gsm and its sin/cos.
    # a[11,...,33] - matrix converts geo to gsm.
    # cgst/sgst - cos/sin of gst.
    # ds3.
    # ba(6).
    global st0,ct0,sl0,cl0,ctcl,stcl,ctsl,stsl,sfi,cfi,sps,cps, \
        shi,chi,hi,psi,xmut,a11,a21,a31,a12,a22,a32,a13,a23,a33,ds3,cgst,sgst,ba

    # The common block /geopack2/ contains coefficients of the IGRF field model, calculated
    # for a given year and day from their standard epoch values. the array rec contains
    # coefficients used in the recursion relations for legendre associate polynomials.
    # common /geopack2/ g(105),h(105),rec(105)
    global g,h,rec


    # Compute the m,n related coefficients (following the notation in Davis 2004):
    # 1. The Schmidt quasi-normalization: S_n,m, which normalizes the associated Legendre polynomials P_n^m
    #    to the Guassian normalized associated Legendre polynomials P^n,m.
    #
    #    Since g_n^m * P_n^m should be constant, mutiplying S_n,m to g_n^m is equivalently converting P_n^m to P^n,m.
    #    The benefit of doing so is that P^n,m follows simple recursive form, c.f. Eq (16 a-c) and Eq(17 a-b).
    # 2. rec[mn], which used in the recursion relation.
    k = 14
    nmn = np.int32((k+1)*k/2)

    # rec[mn].
    rec = np.empty(nmn,dtype=float)
    mn = 0
    for n in range(k):                  # K^1,m = 0, Eq (17a), automatically done.
        n2 = 2*n+1
        n2 = n2*(n2-2)
        for m in range(n+1):
            rec[mn] = (n-m)*(n+m)/n2    # K^n,m = (n-m)(n+m)/(2n+1)(2n-1), Eq (17b)
            mn += 1

    # coefficients for a given time, g_n^m(t), h_n^m(t)
    g,h = load_igrf(ut)

    # now multiply them by schmidt normalization factors:
    s = 1.                      # S_0,0 = 1, Eq (18a)
    mn = 0
    for n in range(1,k):
        mn += 1                 # skip m=0.
        s *= (2*n-1)/n
        g[mn] *= s              # S_n,0 = S_n-1,0 * (2n-1)/n, Eq (18b)
        h[mn] *= s
        p = s
        for m in range(1,n+1):
            if m == 1: aa = 2   # aa = delta_m,1
            else: aa = 1
            p *= np.sqrt(aa*(n-m+1)/(n+m))
            mn += 1
            g[mn] *= p          # S_n,m = S_n,m-1 * sqrt(aa(n-m+1)/(n+m)), Eq (18c)
            h[mn] *= p          # now g/h are actually g^n,m, Eq (14 a-b)

    g10=-g[1]
    g11= g[2]
    h11= h[2]

    # Now calculate the components of the unit vector ezmag in geo coord.system:
    # sin(teta0)*cos(lambda0), sin(teta0)*sin(lambda0), and cos(teta0)
    #       st0 * cl0                st0 * sl0                ct0
    sq=g11**2+h11**2
    sqq=np.sqrt(sq)
    sqr=np.sqrt(g10**2+sq)
    sl0=-h11/sqq
    cl0=-g11/sqq
    st0=sqq/sqr
    ct0=g10/sqr
    stcl=st0*cl0
    stsl=st0*sl0
    ctsl=ct0*sl0
    ctcl=ct0*cl0


    gst,slong,srasn,sdec,obliq = sun(ut)
    # s1,s2, and s3 are the components of the unit vector exgsm=exgse in the
    # system gei pointing from the earth's center to the sun:
    s1=np.cos(srasn)*np.cos(sdec)
    s2=np.sin(srasn)*np.cos(sdec)
    s3=np.sin(sdec)
    cgst=np.cos(gst)
    sgst=np.sin(gst)

    # dip1, dip2, and dip3 are the components of the unit vector ezsm=ezmag in the system gei:
    dip1=stcl*cgst-stsl*sgst
    dip2=stcl*sgst+stsl*cgst
    dip3=ct0

    # Now calculate the components of the unit vector eygsm in the system gei
    # by taking the vector product d x s and normalizing it to unit length:
    y1=dip2*s3-dip3*s2
    y2=dip3*s1-dip1*s3
    y3=dip1*s2-dip2*s1
    y=np.sqrt(y1*y1+y2*y2+y3*y3)
    y1=y1/y
    y2=y2/y
    y3=y3/y

    # Then in the gei system the unit vector z = ezgsm = exgsm x eygsm = s x y has the components:
    z1=s2*y3-s3*y2
    z2=s3*y1-s1*y3
    z3=s1*y2-s2*y1

    # The vector ezgse (here dz) in gei has the components (0,-sin(delta),
    # cos(delta)) = (0.,-0.397823,0.917462); here delta = 23.44214 deg for
    # The epoch 1978 (see the book by gurevich or other astronomical handbooks).
    # Here the most accurate time-dependent formula is used:
    dz1=0.
    dz2=-np.sin(obliq)
    dz3= np.cos(obliq)

    # then the unit vector eygse in gei system is the vector product dz x s :
    dy1=dz2*s3-dz3*s2
    dy2=dz3*s1-dz1*s3
    dy3=dz1*s2-dz2*s1

    # The elements of the matrix gse to gsm are the scalar products:
    # chi=em22=(eygsm,eygse), shi=em23=(eygsm,ezgse), em32=(ezgsm,eygse)=-em23,
    # and em33=(ezgsm,ezgse)=em22
    chi=y1*dy1+y2*dy2+y3*dy3
    shi=y1*dz1+y2*dz2+y3*dz3
    hi=np.arcsin(shi)

    # Tilt angle: psi=arcsin(dip,exgsm)
    sps=dip1*s1+dip2*s2+dip3*s3
    cps=np.sqrt(1.-sps**2)
    psi=np.arcsin(sps)

    # The elements of the matrix mag to sm are the scalar products:
    # cfi=gm22=(eysm,eymag), sfi=gm23=(eysm,exmag); They can be derived as follows:
    # In geo the vectors exmag and eymag have the components (ct0*cl0,ct0*sl0,-st0) and (-sl0,cl0,0), respectively.

    # Hence, in gei the components are:
    #     exmag:    ct0*cl0*cos(gst)-ct0*sl0*sin(gst)
    #               ct0*cl0*sin(gst)+ct0*sl0*cos(gst)
    #               -st0
    #     eymag:    -sl0*cos(gst)-cl0*sin(gst)
    #               -sl0*sin(gst)+cl0*cos(gst)
    #                0
    # The components of eysm in gei were found above as y1, y2, and y3;
    # Now we only have to combine the quantities into scalar products:
    exmagx=ct0*(cl0*cgst-sl0*sgst)
    exmagy=ct0*(cl0*sgst+sl0*cgst)
    exmagz=-st0
    eymagx=-(sl0*cgst+cl0*sgst)
    eymagy=-(sl0*sgst-cl0*cgst)
    cfi=y1*eymagx+y2*eymagy
    sfi=y1*exmagx+y2*exmagy+y3*exmagz

    xmut=(np.arctan2(sfi,cfi)+3.1415926536)*3.8197186342

    # The elements of the matrix geo to gsm are the scalar products:
    # a11=(exgeo,exgsm), a12=(eygeo,exgsm), a13=(ezgeo,exgsm),
    # a21=(exgeo,eygsm), a22=(eygeo,eygsm), a23=(ezgeo,eygsm),
    # a31=(exgeo,ezgsm), a32=(eygeo,ezgsm), a33=(ezgeo,ezgsm),

    # All the unit vectors in brackets are already defined in gei:
    # exgeo=(cgst,sgst,0), eygeo=(-sgst,cgst,0), ezgeo=(0,0,1)
    # exgsm=(s1,s2,s3),  eygsm=(y1,y2,y3),   ezgsm=(z1,z2,z3)
    # and therefore:
    a11=s1*cgst+s2*sgst
    a12=-s1*sgst+s2*cgst
    a13=s3
    a21=y1*cgst+y2*sgst
    a22=-y1*sgst+y2*cgst
    a23=y3
    a31=z1*cgst+z2*sgst
    a32=-z1*sgst+z2*cgst
    a33=z3

    return psi

def sun(ut):
    """
    Calculates four quantities necessary for coordinate transformations
    which depend on sun position (and, hence, on universal time and season)
    Based on http://aa.usno.navy.mil/faq/docs/SunApprox.php and http://aa.usno.navy.mil/faq/docs/GAST.php

    :param ut: ut sec, can be array.
    :return: gst,slong,srasn,sdec. gst - greenwich mean sidereal time, slong - longitude along ecliptic
        srasn - right ascension,  sdec - declination of the sun (radians)
        obliq - mean oblique of the ecliptic (radian)

    Python version by Sheng Tian
    """

    twopi = 2*np.pi
    jd2000 = 2451545.0
    # convert to Julian date.
    t0_jd = 2440587.5       # in day, 0 of Julian day.
    secofday1 = 1./86400    # 1/sec of day.
    t_jd = ut*secofday1+t0_jd

    # d = mjd - mj2000.
    d = t_jd-jd2000
    d = np.squeeze(d)

    # mean obliquity of the ecliptic, e.
    # e = 23.439 - 0.00000036*d     ; in degree.
    e = 0.4090877233749509 - 6.2831853e-9*d

    # mean anomaly of the Sun, g.
    # g = 357.529 + 0.98560028*d    ; in degree.
    g = 6.2400582213628066 + 0.0172019699945780*d
    g = np.mod(g, twopi)

    # mean longitude of the Sun, q.
    # q = 280.459 + 0.98564736*d   ; in degree.
    q = 4.8949329668507771 + 0.0172027916955899*d
    q = np.mod(q, twopi)

    # geocentric apparent ecliptic longitude, l.
    # l = q + 1.915 sin g + 0.020 sin 2g    ; in degree.
    l = q + 0.0334230551756914*np.sin(g) + 0.0003490658503989*np.sin(2*g)

    # vl - q, mean longitude of the sun.
    # vl = np.mod(279.696678+0.9856473354*dj,360.)/rad
    # q = np.mod(4.881627937990388+0.01720279126623886*dj, twopi)

    # g, mean anomaly of the sun.
    # g = np.mod(358.475845+0.985600267*dj,360.)/rad
    # g = np.mod(6.256583784118852+0.017201969767685215*dj, twopi)

    # slong - l, geocentric apparent ecliptic longitude.
    # slong = (vl + (1.91946-0.004789*t)*sin(g) + 0.020094*sin(2*g))/rad
    # l = q+(0.03350089686033036-2.2884002156881157e-09*dj)*np.sin(g)+0.0003507064598957406*np.sin(2*g)
    # l = np.mod(l, twopi)

    # obliq - e, mean obliquity of the ecliptic.
    # obliq = (23.45229-0.0130125*t)/rad
    # e = 0.40931967763254096-6.217959450123535e-09*dj

    # sin(d) = sin(e) * sin(L)
    sind = np.sin(e)*np.sin(l)
    sdec = np.arcsin(sind)

    # tan(RA) = cos(e)*sin(L)/cos(L)
    srasn = np.arctan2(np.cos(e)*np.sin(l), np.cos(l))
    srasn = np.mod(srasn, twopi)


    # http://aa.usno.navy.mil/faq/docs/GAST.php
    # gst - gmst, greenwich mean sidereal time.
    # gst = np.mod(279.690983+.9856473354*dj+360.*fday+180.,360.)/rad
    # gst = np.mod(4.881528541489487+0.01720279126623886*dj+twopi*fday+np.pi, twopi)
    # gmst = 18.697374558 + 24.06570982441908*d  # in hour
    gmst = 4.894961212735792 + 6.30038809898489*d  # in rad
    gmst = np.mod(gmst, twopi)

    return gmst,l,srasn,sdec,e



def geomag(p1,p2,p3, j):
    """
    Converts geographic (geo) to dipole (mag) coordinates or vica versa.
                   j>0                       j<0
    input:  j,xgeo,ygeo,zgeo           j,xmag,ymag,zmag
    output:    xmag,ymag,zmag           xgeo,ygeo,zgeo

    :param p1,p2,p3: input position
    :param j: flag
    :return: output position
    """

    #  Attention: subroutine  recalc  must be invoked before geomag in two cases:
    #     /a/  before the first transformation of coordinates
    #     /b/  if the values of time have been changed

    # common /geopack1/ st0,ct0,sl0,cl0,ctcl,stcl,ctsl,stsl,ab(19),bb(8)
    global st0,ct0, sl0,sl0, ctcl,stcl, ctsl,stsl

    if j > 0:
        xgeo,ygeo,zgeo = [p1,p2,p3]
        xmag = xgeo*ctcl+ygeo*ctsl-zgeo*st0
        ymag = ygeo*cl0-xgeo*sl0
        zmag = xgeo*stcl+ygeo*stsl+zgeo*ct0
        return xmag,ymag,zmag
    else:
        xmag,ymag,zmag = [p1,p2,p3]
        xgeo = xmag*ctcl-ymag*sl0+zmag*stcl
        ygeo = xmag*ctsl+ymag*cl0+zmag*stsl
        zgeo = zmag*ct0-xmag*st0
        return xgeo,ygeo,zgeo

def geigeo(p1,p2,p3, j):
    """
    Converts equatorial inertial (gei) to geographical (geo) coords or vica versa.
                   j>0                       j<0
    input:  j,xgei,ygei,zgei           j,xgeo,ygeo,zgeo
    output:    xgeo,ygeo,zgeo           xgei,ygei,zgei

    :param p1,p2,p3: input position
    :param j: flag
    :return: output position
    """

    #  Attention: subroutine  recalc  must be invoked before geomag in two cases:
    #     /a/  before the first transformation of coordinates
    #     /b/  if the values of time have been changed

    # common /geopack1/ a(27),cgst,sgst,b(6)
    global cgst,sgst

    if j > 0:
        xgei,ygei,zgei = [p1,p2,p3]
        xgeo = xgei*cgst+ygei*sgst
        ygeo = ygei*cgst-xgei*sgst
        zgeo = zgei
        return xgeo,ygeo,zgeo
    else:
        xgeo,ygeo,zgeo = [p1,p2,p3]
        xgei = xgeo*cgst-ygeo*sgst
        ygei = ygeo*cgst+xgeo*sgst
        zgei = zgeo
        return xgei,ygei,zgei

def magsm(p1,p2,p3, j):
    """
    Converts dipole (mag) to solar magnetic (sm) coordinates or vica versa
                   j>0                       j<0
    input:  j,xmag,ymag,zmag           j,xsm, ysm, zsm
    output:    xsm, ysm, zsm           xmag,ymag,zmag

    :param p1,p2,p3: input position
    :param j: flag
    :return: output position
    """

    #  Attention: subroutine  recalc  must be invoked before geomag in two cases:
    #     /a/  before the first transformation of coordinates
    #     /b/  if the values of time have been changed

    # common /geopack1/ a(8),sfi,cfi,b(7),ab(10),ba(8)
    global sfi,cfi

    if j > 0:
        xmag,ymag,zmag = [p1,p2,p3]
        xsm = xmag*cfi-ymag*sfi
        ysm = xmag*sfi+ymag*cfi
        zsm = zmag
        return xsm,ysm,zsm
    else:
        xsm,ysm,zsm = [p1,p2,p3]
        xmag = xsm*cfi+ysm*sfi
        ymag = ysm*cfi-xsm*sfi
        zmag = zsm
        return xmag,ymag,zmag

def gsmgse(p1,p2,p3, j):
    """
    converts geocentric solar magnetospheric (gsm) coords to solar ecliptic (gse) ones or vica versa.
                   j>0                       j<0
    input:  j,xgsm,ygsm,zgsm           j,xgse,ygse,zgse
    output:    xgse,ygse,zgse           xgsm,ygsm,zgsm

    :param p1,p2,p3: input position
    :param j: flag
    :return: output position
    """

    # common /geopack1/ a(12),shi,chi,ab(13),ba(8)
    global shi,chi

    if j > 0:
        xgsm,ygsm,zgsm = [p1,p2,p3]
        xgse = xgsm
        ygse = ygsm*chi-zgsm*shi
        zgse = ygsm*shi+zgsm*chi
        return xgse,ygse,zgse
    else:
        xgse,ygse,zgse = [p1,p2,p3]
        xgsm = xgse
        ygsm = ygse*chi+zgse*shi
        zgsm = zgse*chi-ygse*shi
        return xgsm,ygsm,zgsm

def smgsm(p1,p2,p3, j):
    """
    Converts solar magnetic (sm) to geocentric solar magnetospheric (gsm) coordinates or vica versa.
                   j>0                       j<0
    input:  j,xsm, ysm, zsm           j,xgsm,ygsm,zgsm
    output:    xgsm,ygsm,zgsm           xsm, ysm, zsm

    :param p1,p2,p3: input position
    :param j: flag
    :return: output position
    """

    #  Attention: subroutine  recalc  must be invoked before geomag in two cases:
    #     /a/  before the first transformation of coordinates
    #     /b/  if the values of time have been changed

    # common /geopack1/ a(10),sps,cps,b(15),ab(8)
    global sps,cps

    if j > 0:
        xsm,ysm,zsm = [p1,p2,p3]
        xgsm = xsm*cps+zsm*sps
        ygsm = ysm
        zgsm = zsm*cps-xsm*sps
        return xgsm,ygsm,zgsm
    else:
        xgsm,ygsm,zgsm = [p1,p2,p3]
        xsm = xgsm*cps-zgsm*sps
        ysm = ygsm
        zsm = xgsm*sps+zgsm*cps
        return xsm,ysm,zsm

def geogsm(p1,p2,p3, j):
    """
    Converts geographic (geo) to geocentric solar magnetospheric (gsm) coordinates or vica versa.
                   j>0                       j<0
    input:  j,xgeo,ygeo,zgeo           j,xgsm,ygsm,zgsm
    output:    xgsm,ygsm,zgsm           xgeo,ygeo,zgeo

    :param p1,p2,p3: input position
    :param j: flag
    :return: output position
    """

    #  Attention: subroutine  recalc  must be invoked before geomag in two cases:
    #     /a/  before the first transformation of coordinates
    #     /b/  if the values of time have been changed

    # common /geopack1/aa(17),a11,a21,a31,a12,a22,a32,a13,a23,a33,d,b(8)
    global a11,a21,a31,a12,a22,a32,a13,a23,a33

    if j > 0:
        xgeo,ygeo,zgeo = [p1,p2,p3]
        xgsm = a11*xgeo+a12*ygeo+a13*zgeo
        ygsm = a21*xgeo+a22*ygeo+a23*zgeo
        zgsm = a31*xgeo+a32*ygeo+a33*zgeo
        return xgsm,ygsm,zgsm
    else:
        xgsm,ygsm,zgsm = [p1,p2,p3]
        xgeo = a11*xgsm+a21*ygsm+a31*zgsm
        ygeo = a12*xgsm+a22*ygsm+a32*zgsm
        zgeo = a13*xgsm+a23*ygsm+a33*zgsm
        return xgeo,ygeo,zgeo



def sphcar(p1,p2,p3, j):
    """
    Converts spherical coords into cartesian ones and vica versa (theta and phi in radians).
                  j>0            j<0
    input:   j,r,theta,phi     j,x,y,z
    output:      x,y,z        r,theta,phi

    :param r,theta,phi:
    :param x,y,z:
    :param j:
    :return:

    Note: at the poles (x=0 and y=0) we assume phi=0 (when converting from cartesian to spherical coords, i.e., for j<0)
    Last mofification: April 1, 2003 (only some notation changes and more comments added)
    Author:  N.A. Tsyganenko
    """

    if j > 0:
        r,theta,phi = [p1,p2,p3]
        sq=r*np.sin(theta)
        x=sq*np.cos(phi)
        y=sq*np.sin(phi)
        z= r*np.cos(theta)
        return x,y,z
    else:
        x,y,z = [p1,p2,p3]
        sq=x**2+y**2
        r=np.sqrt(sq+z**2)
        if sq != 0:
            sq=np.sqrt(sq)
            phi=np.arctan2(y,x)
            theta=np.arctan2(sq,z)
            if phi < 0: phi += 2*np.pi
        else:
            phi=0.
            if z < 0: theta = np.pi
            else: theta = 0
        return r,theta,phi

def bspcar(theta,phi,br,btheta,bphi):
    """
    Calculates cartesian field components from spherical ones.

    :param theta,phi: spherical angles of the point in radians
    :param br,btheta,bphi: spherical components of the field
    :return: bx,by,bz. cartesian components of the field

    # Last mofification:  April 1, 2003 (only some notation changes and more comments added)
    # Author:  N.A. Tsyganenko
    """

    s= np.sin(theta)
    c= np.cos(theta)
    sf=np.sin(phi)
    cf=np.cos(phi)
    be=br*s+btheta*c

    bx=be*cf-bphi*sf
    by=be*sf+bphi*cf
    bz=br*c-btheta*s

    return bx,by,bz

def bcarsp(x,y,z,bx,by,bz):
    """
    Calculates spherical field components from those in cartesian system

    :param x,y,z: cartesian components of the position vector
    :param bx,by,bz: cartesian components of the field vector
    :return: br,btheta,bphi. spherical components of the field vector

    Last mofification:  April 1, 2003 (only some notation changes and more comments added)
    Author:  N.A. Tsyganenko
    """

    rho2=x**2+y**2
    r=np.sqrt(rho2+z**2)
    rho=np.sqrt(rho2)

    if rho != 0:
        cphi=x/rho
        sphi=y/rho
    else:
        cphi=1.
        sphi=0.

    ct=z/r
    st=rho/r

    br=(x*bx+y*by+z*bz)/r
    btheta=(bx*cphi+by*sphi)*ct-bz*st
    bphi=by*cphi-bx*sphi

    return br,btheta,bphi


def call_external_model(exname, par, ps, x,y,z):
    if exname == 't89':
        return t89.t89(par, ps, x, y, z)
    elif exname == 't96':
        return t96.t96(par, ps, x, y, z)
    elif exname == 't01':
        return t01.t01(par, ps, x, y, z)
    elif exname == 't04':
        return t04.t04(par, ps, x, y, z)
    else:
        raise ValueError

def call_internal_model(inname, x,y,z):
    if inname == 'dipole':
        return dip(x,y,z)
    elif inname == 'igrf':
        return igrf_gsm(x,y,z)
    else:
        raise ValueError

def rhand(x,y,z,parmod,exname,inname):
    """
    Calculates the components of the right hand side vector in the geomagnetic field
    line equation  (a subsidiary subroutine for the subroutine step)

    :param x,y,z:
    :param parmod:
    :param exname: name of the subroutine for the external field.
    :param inname: name of the subroutine for the internal field.

    Last mofification:  March 31, 2003
    Author:  N.A. Tsyganenko
    :return: r1,r2,r3.
    """
    #  common /geopack1/ a(15),psi,aa(10),ds3,bb(8)
    global a, psi, aa, ds3, bb

    bxgsm,bygsm,bzgsm = call_external_model(exname, parmod, psi, x,y,z)
    hxgsm,hygsm,hzgsm = call_internal_model(inname, x,y,z)

    bx=bxgsm+hxgsm
    by=bygsm+hygsm
    bz=bzgsm+hzgsm
    b=ds3/np.sqrt(bx**2+by**2+bz**2)

    r1=bx*b
    r2=by*b
    r3=bz*b

    return r1,r2,r3

def step(x,y,z, ds,errin,parmod,exname,inname):
    """
    Re-calculates {x,y,z}, making a step along a field line.
    model version, the array parmod contains input parameters for that model

    :param x,y,z: the input position
    :param ds: the step size
    :param errin: the permissible error value
    :param parmod: contains input parameters for that model
    :param exname: name of the subroutine for the external field.
    :param inname: name of the subroutine for the internal field.
    :return: x,y,z. The output position

    Last mofification:  March 31, 2003
    Author:  N.A. Tsyganenko
    """

    # common /geopack1/ a(26),ds3,b(8)
    global a, ds3, b

    if errin <=0: raise ValueError
    errcur = errin
    i = 0
    maxloop = 100

    while (errcur >= errin) & (i < maxloop):
        ds3=-ds/3.
        r11,r12,r13 = rhand(x,y,z,parmod,exname,inname)
        r21,r22,r23 = rhand(x+r11,y+r12,z+r13,parmod,exname,inname)
        r31,r32,r33 = rhand(x+.5*(r11+r21),y+.5*(r12+r22),z+.5*(r13+r23),parmod,exname,inname)
        r41,r42,r43 = rhand(x+.375*(r11+3.*r31),y+.375*(r12+3.*r32),z+.375*(r13+3.*r33),parmod,exname,inname)
        r51,r52,r53 = rhand(x+1.5*(r11-3.*r31+4.*r41),y+1.5*(r12-3.*r32+4.*r42),z+1.5*(r13-3.*r33+4.*r43),parmod,exname,inname)
        errcur=np.abs(r11-4.5*r31+4.*r41-.5*r51)+np.abs(r12-4.5*r32+4.*r42-.5*r52)+np.abs(r13-4.5*r33+4.*r43-.5*r53)
        if errcur < errin: break

        ds *= 0.5
        i += 1
    else:
        print('reached maximum loop ...')
        return

    x += 0.5*(r11+4.*r41+r51)
    y += 0.5*(r12+4.*r42+r52)
    z += 0.5*(r13+4.*r43+r53)
    if (errcur < (errin*0.04)) & (np.abs(ds) < 1.33):
        ds *= 1.5

    return x,y,z

def trace(xi,yi,zi,dir,rlim=10,r0=1,parmod=2,exname='t89',inname='igrf'):
    """
    Traces a field line from an arbitrary point of space to the earth's surface or
    to a model limiting boundary.

    The highest order of spherical harmonics in the main field expansion used
    in the mapping is calculated automatically. if inname=igrf_gsm, then an IGRF model
    field will be used, and if inname=dip, a pure dipole field will be used.

    In any case, before calling trace, one should invoke recalc, to calculate correct
    values of the IGRF coefficients and all quantities needed for transformations
    between coordinate systems involved in this calculations.

    Alternatively, the subroutine recalc can be invoked with the desired values of
    ut sec (to specify the dipole moment), while the values of the dipole
    tilt angle psi (in radians) and its sine (sps) and cosine (cps) can be explicitly
    specified and forwarded to the common block geopack1 (11th, 12th, and 16th elements, resp.)

    :param xi,yi,zi: gsm coords of initial point (in earth radii, 1 re = 6371.2 km)
    :param dir: sign of the tracing direction: if
        dir=1.0 then we move antiparallel to the field vector (e.g. from northern to southern conjugate point), and if
        dir=-1.0 then the tracing goes in the opposite direction.
    :param rlim: upper limit of the geocentric distance, where the tracing is terminated.
    :param r0: radius of a sphere (in re) for which the field line endpoint coordinates xf,yf,zf should be calculated.
    :param parmod: the iopt for T89, or for the other Txx models it is
        a 10-element array containing model parameters, needed for a unique specification of the external field.
        The concrete meaning of the components of parmod depends on a specific version of the external field model.
    :param exname: name of the subroutine for the external field.
    :param inname: name of the subroutine for the internal field.
    :return:
        xf,yf,zf. GSM coords of the last calculated point of a field line
        xx,yy,zz. Arrays containing coords of field line points. Here their maximal length was assumed equal to 999.
        l. Actual number of the calculated field line points. if l exceeds 999, tracing terminates, and a warning is displayed.

    Last mofification:  March 31, 2003
    Author:  N.A. Tsyganenko
    """

    # common /geopack1/ aa(26),dd,bb(8)
    global aa, dd, bb, ds3

    err, l, ds, x,y,z, dd, al = 0.001, 0, 0.5*dir, xi,yi,zi, dir, 0.
    xx = np.array([x])
    yy = np.array([y])
    zz = np.array([z])
    maxloop = 1000

    # Here we call RHAND just to find out the sign of the radial component of the field
    # vector, and to determine the initial direction of the tracing (i.e., either away
    # or towards Earth):
    ds3 = -ds/3.
    r1,r2,r3 = rhand(x,y,z,parmod,exname,inname)

    # |ad|=0.01 and its sign follows the rule:
    # (1) if dir=1 (tracing antiparallel to B vector) then the sign of ad is the same as of Br
    # (2) if dir=-1 (tracing parallel to B vector) then the sign of ad is opposite to that of Br
    #     ad is defined in order to initialize the value of rr (radial distance at previous step):
    ad=0.01
    if (x*r1+y*r2+z*r3) < 0: ad=-0.01

    rr=np.sqrt(x**2+y**2+z**2)+ad

    while l < maxloop:
        ryz=y**2+z**2
        r2=x**2+ryz
        r=np.sqrt(r2)

        # check if the line hit the outer tracing boundary; if yes, then terminate the tracing
        if (r > rlim) | (ryz > 1600) | (x>20): break

        # check whether or not the inner tracing boundary was crossed from outside,
        # if yes, then calculate the footpoint position by interpolation
        if (r < r0) & (rr > r):
            # find the footpoint position by interpolating between the current and previous field line points:
            r1=(r0-r)/(rr-r)
            x=x-(x-xr)*r1
            y=y-(y-yr)*r1
            z=z-(z-zr)*r1
            break

        # check if (i) we are moving outward, or (ii) we are still sufficiently
        # far from Earth (beyond R=5Re); if yes, proceed further:
        if (r >= rr) | (r > 5):
            pass
        # now we moved closer inward (between R=3 and R=5); go to 3 and begin logging
        # previous values of X,Y,Z, to be used in the interpolation (after having
        # crossed the inner tracing boundary):
        else:
            if r >= 3:
                # We entered inside the sphere R=3: to avoid too large steps (and hence inaccurate
                # interpolated position of the footpoint), enforce the progressively smaller
                # stepsize values as we approach the inner boundary R=R0:
                ds = dir
            else:
                fc = 0.2
                if (r-r0) < 0.05: fc = 0.05
                al = fc*(r-r0+0.2)
                ds = dir*al
            xr,yr,zr = [x,y,z]
        rr=r
        x,y,z = step(x,y,z,ds,err,parmod,exname,inname)
        np.append(xx,x)
        np.append(yy,y)
        np.append(zz,z)
        l += 1

    return x,y,z



def shuetal_mgnp(xn_pd,vel,bzimf,xgsm,ygsm,zgsm):
    """
    For any point of space with coordinates (xgsm,ygsm,zgsm) and specified conditions
    in the incoming solar wind, this subroutine:
        (1) determines if the point (xgsm,ygsm,zgsm) lies inside or outside the
            model magnetopause of Shue et al. (jgr-a, v.103, p. 17691, 1998).
        (2) calculates the gsm position of a point {xmgnp,ymgnp,zmgnp}, lying at the model
            magnetopause and asymptotically tending to the nearest boundary point with
            respect to the observation point {xgsm,ygsm,zgsm}, as it approaches the magnetopause.

    :param xn_pd: either solar wind proton number density (per c.c.) (if vel>0)
        or the solar wind ram pressure in nanopascals   (if vel<0)
    :param vel: either solar wind velocity (km/sec)
        or any negative number, which indicates that xn_pd stands
        for the solar wind pressure, rather than for the density
    :param bzimf: imf bz in nanoteslas
    :param xgsm,ygsm,zgsm: gsm position of the observation point in earth radii
    :return: xmgnp,ymgnp,zmgnp. GSM position of the boundary point.
        dist. Distance (in re) between the observation point (xgsm,ygsm,zgsm) and the model magnetopause
        id. position flag: id=+1 (-1) means that the observation point lies inside (outside) of the model magnetopause, respectively.

    Last mofification:  April 4, 2003
    Author:  N.A. Tsyganenko
    """

    # pd is the solar wind dynamic pressure (in npa)
    if vel < 0: pd = xn_pd
    else: pd = 1.94e-6*xn_pd*vel**2

    # Define the angle phi, measured duskward from the noon-midnight meridian plane;
    # if the observation point lies on the x axis, the angle phi cannot be uniquely
    # defined, and we set it at zero:
    if (ygsm != 0) | (zgsm != 0): phi = np.arctan2(ygsm,zgsm)
    else: phi = 0

    # First, find out if the observation point lies inside the Shue et al bdry
    # and set the value of the id flag:
    id = -1
    r0 = (10.22+1.29*np.tanh(0.184*(bzimf+8.14)))*pd**(-0.15151515)
    alpha = (0.58-0.007*bzimf)*(1.+0.024*np.log(pd))
    r = np.sqrt(xgsm**2+ygsm**2+zgsm**2)
    rm = r0*(2./(1.+xgsm/r))**alpha
    if r < rm: id = 1

    #  Now, find the corresponding t96 magnetopause position, to be used as
    #  a starting approximation in the search of a corresponding Shue et al.
    #  boundary point:
    xmt96,ymt96,zmt96,dist,id96 = t96_mgnp(pd,-1.,xgsm,ygsm,zgsm)

    rho2 = ymt96**2+zmt96**2
    r = np.sqrt(rho2+xmt96**2)
    st = np.sqrt(rho2)/r
    ct = xmt96/r

    #  Now, use newton's iterative method to find the nearest point at the Shue et al.'s boundary:
    nit = 0
    while nit < 1000:
        nit += 1

        t = np.arctan2(st,ct)
        rm = r0*(2./(1.+ct))**alpha

        f = r-rm
        gradf_r = 1.
        gradf_t = -alpha/r*rm*st/(1.+ct)
        gradf = np.sqrt(gradf_r**2+gradf_t**2)

        dr = -f/gradf**2
        dt =  dr/r*gradf_t

        r = r+dr
        t = t+dt
        st = np.sin(t)
        ct = np.cos(t)

        ds = np.sqrt(dr**2+(r*dt)**2)
        if ds <= 1.e-4: break
    else: print(' boundary point could not be found; iterations do not converge')


    xmgnp = r*np.cos(t)
    rho =   r*np.sin(t)
    ymgnp = rho*np.sin(phi)
    zmgnp = rho*np.cos(phi)

    dist = np.sqrt((xgsm-xmgnp)**2+(ygsm-ymgnp)**2+(zgsm-zmgnp)**2)

    return xmgnp,ymgnp,zmgnp, dist, id

def t96_mgnp(xn_pd,vel,xgsm,ygsm,zgsm):
    """
    For any point of space with given coordinates (xgsm,ygsm,zgsm), this subroutine defines
    the position of a point (xmgnp,ymgnp,zmgnp) at the T96 model magnetopause, having the
    same value of the ellipsoidal tau-coordinate, and the distance between them. This is
    not the shortest distance d_min to the boundary, but dist asymptotically tends to d_min,
    as the observation point gets closer to the magnetopause.

    The pressure-dependent magnetopause is that used in the t96_01 model
    (Tsyganenko, jgr, v.100, p.5599, 1995; esa sp-389, p.181, oct. 1996)

    :param xn_pd: either solar wind proton number density (per c.c.) (if vel>0)
        or the solar wind ram pressure in nanopascals   (if vel<0)
    :param vel: either solar wind velocity (km/sec)
        or any negative number, which indicates that xn_pd stands
        for the solar wind pressure, rather than for the density
    :param xgsm,ygsm,zgsm: gsm position of the observation point in earth radii
    :return: xmgnp,ymgnp,zmgnp. GSM position of the boundary point.
        dist. Distance (in re) between the observation point (xgsm,ygsm,zgsm) and the model magnetopause
        id. position flag: id=+1 (-1) means that the observation point lies inside (outside) of the model magnetopause, respectively.

    Last mofification:  April 3, 2003
    Author:  N.A. Tsyganenko
    """
    # Define solar wind dynamic pressure (nanopascals, assuming 4% of alpha-particles),
    # if not explicitly specified in the input:
    if vel < 0: pd = xn_pd
    else: pd = 1.94e-6*xn_pd*vel**2

    # ratio of pd to the average pressure, assumed equal to 2 npa:
    # The power index 0.14 in the scaling factor is the best-fit value
    # obtained from data and used in the t96_01 version
    # Values of the magnetopause parameters for  pd = 2 npa:
    rat = pd/2.0
    rat16 = rat**0.14

    a0, s00, x00 = [70.,1.08,5.48]

    # Values of the magnetopause parameters, scaled by the actual pressure:
    # xm is the x-coordinate of the "seam" between the ellipsoid and the cylinder
    # For details on the ellipsoidal coordinates, see the paper:
    # N.A. Tsyganenko, Solution of chapman-ferraro problem for an ellipsoidal magnetopause, planet.space sci., v.37, p.1037, 1989).
    a  = a0/rat16
    s0 = s00
    x0 = x00/rat16
    xm = x0-a


    if (ygsm != 0) | (zgsm != 0): phi = np.arctan2(ygsm,zgsm)
    else: phi = 0

    rho = np.sqrt(ygsm**2+zgsm**2)
    if xgsm < xm:
        xmgnp = xgsm
        rhomgnp = a*np.sqrt(s0**2-1)
        ymgnp = rhomgnp*np.sin(phi)
        zmgnp = rhomgnp*np.cos(phi)
        dist = np.sqrt((xgsm-xmgnp)**2+(ygsm-ymgnp)**2+(zgsm-zmgnp)**2)
        if rhomgnp >  rho: id =  1
        else: id = -1
        return xmgnp,ymgnp,zmgnp, dist, id

    xksi = (xgsm-x0)/a+1.
    xdzt = rho/a
    sq1 = np.sqrt((1.+xksi)**2+xdzt**2)
    sq2 = np.sqrt((1.-xksi)**2+xdzt**2)
    sigma = 0.5*(sq1+sq2)
    tau   = 0.5*(sq1-sq2)

    # Now calculate (x,y,z) for the closest point at the magnetopause
    xmgnp = x0-a*(1.-s0*tau)
    arg = (s0**2-1.)*(1.-tau**2)
    if arg < 0: arg = 0
    rhomgnp = a*np.sqrt(arg)
    ymgnp = rhomgnp*np.sin(phi)
    zmgnp = rhomgnp*np.cos(phi)

    # Now calculate the distance between the points {xgsm,ygsm,zgsm} and {xmgnp,ymgnp,zmgnp}:
    # (in general, this is not the shortest distance d_min, but dist asymptotically tends
    # to d_min, as we are getting closer to the magnetopause):
    dist = np.sqrt((xgsm-xmgnp)**2+(ygsm-ymgnp)**2+(zgsm-zmgnp)**2)
    if sigma > s0: id = -1  # id = -1 means that the point lies outside
    else: id = 1            # id =  1 means that the point lies inside

    return xmgnp,ymgnp,zmgnp, dist, id


init_igrf()
