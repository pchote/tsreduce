#include <math.h>
#include "astro_convert.h"
#include "helpers.h"

// Table of leap seconds for calculating TAI from UTC
// Obtained from ftp://time.nist.gov/pub/leap-seconds.list
// First column is seconds after 1900, second is total leap seconds
double leapseconds[] =
{
    2272060800, 10, // 1 Jan 1972
    2287785600, 11, // 1 Jul 1972
    2303683200, 12, // 1 Jan 1973
    2335219200, 13, // 1 Jan 1974
    2366755200, 14, // 1 Jan 1975
    2398291200, 15, // 1 Jan 1976
    2429913600, 16, // 1 Jan 1977
    2461449600, 17, // 1 Jan 1978
    2492985600, 18, // 1 Jan 1979
    2524521600, 19, // 1 Jan 1980
    2571782400, 20, // 1 Jul 1981
    2603318400, 21, // 1 Jul 1982
    2634854400, 22, // 1 Jul 1983
    2698012800, 23, // 1 Jul 1985
    2776982400, 24, // 1 Jan 1988
    2840140800, 25, // 1 Jan 1990
    2871676800, 26, // 1 Jan 1991
    2918937600, 27, // 1 Jul 1992
    2950473600, 28, // 1 Jul 1993
    2982009600, 29, // 1 Jul 1994
    3029443200, 30, // 1 Jan 1996
    3076704000, 31, // 1 Jul 1997
    3124137600, 32, // 1 Jan 1999
    3345062400, 33, // 1 Jan 2006
    3439756800, 34, // 1 Jan 2009
    3550089600, 35  // 1 Jul 2012
};

#define dmod(a1,a2) ((a1)-(a2)*floor((a1)/(a2)))

double dprema[3][3],dpsi,d1pdro,dsinls,dcosls,dsinep,dcosep,forbel[7];
double sorbel[17],sinlp[4],coslp[4],sinlm,coslm,sigma;
int  ideq;
double dc2pi = 6.2831853071796, dc1 = 1.0, dct0 = 2415020.0, dcjul = 36525.0;
double dcbes = 0.313, dctrop = 365.24219572;

// constants dcfel(i,k) of fast changing elements
// i=1             i=2               i=3
double dcfel[8][3] = {
    {1.7400353e+00, 6.2833195099091e+02, 5.2796e-06},
    {6.2565836e+00, 6.2830194572674e+02,-2.6180e-06},
    {4.7199666e+00, 8.3997091449254e+03,-1.9780e-05},
    {1.9636505e-01, 8.4334662911720e+03,-5.6044e-05},
    {4.1547339e+00, 5.2993466764997e+01, 5.8845e-06},
    {4.6524223e+00, 2.1354275911213e+01, 5.6797e-06},
    {4.2620486e+00, 7.5025342197656e+00, 5.5317e-06},
    {1.4740694e+00, 3.8377331909193e+00, 5.6093e-06}
};

// constants dceps and ccsel(i,k) of slowly changing elements
// i=1           i=2           i=3
double dceps[3] = { 4.093198e-01, -2.271110e-04, -2.860401e-08};
double ccsel[17][3] = {
    {1.675104e-02,-4.179579e-05,-1.260516e-07},
    {2.220221e-01, 2.809917e-02, 1.852532e-05},
    {1.589963e+00, 3.418075e-02, 1.430200e-05},
    {2.994089e+00, 2.590824e-02, 4.155840e-06},
    {8.155457e-01, 2.486352e-02, 6.836840e-06},
    {1.735614e+00, 1.763719e-02, 6.370440e-06},
    {1.968564e+00, 1.524020e-02,-2.517152e-06},
    {1.282417e+00, 8.703393e-03, 2.289292e-05},
    {2.280820e+00, 1.918010e-02, 4.484520e-06},
    {4.833473e-02, 1.641773e-04,-4.654200e-07},
    {5.589232e-02,-3.455092e-04,-7.388560e-07},
    {4.634443e-02,-2.658234e-05, 7.757000e-08},
    {8.997041e-03, 6.329728e-06,-1.939256e-09},
    {2.284178e-02,-9.941590e-05, 6.787400e-08},
    {4.350267e-02,-6.839749e-05,-2.714956e-07},
    {1.348204e-02, 1.091504e-05, 6.903760e-07},
    {3.106570e-02,-1.665665e-04,-1.590188e-07}
};

// constants of the arguments of the short-period perturbations
// by the planets: dcargs(i,k)
// i=1                   i=2
double dcargs[15][2] = {
    {5.0974222e+00,-7.8604195454652e+02},
    {3.9584962e+00,-5.7533848094674e+02},
    {1.6338070e+00,-1.1506769618935e+03},
    {2.5487111e+00,-3.9302097727326e+02},
    {4.9255514e+00,-5.8849265665348e+02},
    {1.3363463e+00,-5.5076098609303e+02},
    {1.6072053e+00,-5.2237501616674e+02},
    {1.3629480e+00,-1.1790629318198e+03},
    {5.5657014e+00,-1.0977134971135e+03},
    {5.0708205e+00,-1.5774000881978e+02},
    {3.9318944e+00, 5.2963464780000e+01},
    {4.8989497e+00, 3.9809289073258e+01},
    {1.3097446e+00, 7.7540959633708e+01},
    {3.5147141e+00, 7.9618578146517e+01},
    {3.5413158e+00,-5.4868336758022e+02}
};

// amplitudes ccamps(n,k) of the short-period perturbations
// n=1          n=2            n=3         n=4         n=5
double ccamps[15][5] = {
    {-2.279594e-5, 1.407414e-5, 8.273188e-6, 1.340565e-5,-2.490817e-7},
    {-3.494537e-5, 2.860401e-7, 1.289448e-7, 1.627237e-5,-1.823138e-7},
    { 6.593466e-7, 1.322572e-5, 9.258695e-6,-4.674248e-7,-3.646275e-7},
    { 1.140767e-5,-2.049792e-5,-4.747930e-6,-2.638763e-6,-1.245408e-7},
    { 9.516893e-6,-2.748894e-6,-1.319381e-6,-4.549908e-6,-1.864821e-7},
    { 7.310990e-6,-1.924710e-6,-8.772849e-7,-3.334143e-6,-1.745256e-7},
    {-2.603449e-6, 7.359472e-6, 3.168357e-6, 1.119056e-6,-1.655307e-7},
    {-3.228859e-6, 1.308997e-7, 1.013137e-7, 2.403899e-6,-3.736225e-7},
    { 3.442177e-7, 2.671323e-6, 1.832858e-6,-2.394688e-7,-3.478444e-7},
    { 8.702406e-6,-8.421214e-6,-1.372341e-6,-1.455234e-6,-4.998479e-8},
    {-1.488378e-6,-1.251789e-5, 5.226868e-7,-2.049301e-7, 0.0e0},
    {-8.043059e-6,-2.991300e-6, 1.473654e-7,-3.154542e-7, 0.0e0},
    { 3.699128e-6,-3.316126e-6, 2.901257e-7, 3.407826e-7, 0.0e0},
    { 2.550120e-6,-1.241123e-6, 9.901116e-8, 2.210482e-7, 0.0e0},
    {-6.351059e-7, 2.341650e-6, 1.061492e-6, 2.878231e-7, 0.0e0}
};

// constants of the secular perturbations in longitude
// ccsec3 and ccsec(n,k)
// n=1           n=2           n=3
double ccsec3 = -7.757020e-08;
double ccsec[4][3] = {
    {1.289600e-06, 5.550147e-01, 2.076942e+00},
    {3.102810e-05, 4.035027e+00, 3.525565e-01},
    {9.124190e-06, 9.990265e-01, 2.622706e+00},
    {9.793240e-07, 5.508259e+00, 1.559103e+01}
};

// sideral rate dcsld in longitude, rate ccsgd in mean anomaly
double dcsld = 1.990987e-07, ccsgd = 1.990969e-07;

// some constants used in the calculation of the lunar contribution
double cckm = 3.122140e-05, ccmld = 2.661699e-06, ccfdi = 2.399485e-07;

// constants dcargm(i,k) of the arguments of the perturbations
// of the motion of the moon
// i=1               i=2
double dcargm[3][2] = {
    {5.1679830e+00, 8.3286911095275e+03},
    {5.4913150e+00,-7.2140632838100e+03},
    {5.9598530e+00, 1.5542754389685e+04}
};

// amplitudes ccampm(n,k) of the perturbations of the moon
// n=1          n=2           n=3           n=4
double ccampm[3][4] = {
    { 1.097594e-01, 2.896773e-07, 5.450474e-02, 1.438491e-07},
    {-2.223581e-02, 5.083103e-08, 1.002548e-02,-2.291823e-08},
    { 1.148966e-02, 5.658888e-08, 8.249439e-03, 4.063015e-08}
};

// ccpamv(k)=a*m*dl/dt (planets), dc1mme=1-mass(earth+moon)
double ccpamv[4] = {8.326827e-11, 1.843484e-11, 1.988712e-12, 1.881276e-12};
double dc1mme=0.99999696;

// ccpam(k)=a*m(planets),ccim=inclination(moon),dc1mme=1-mass(emb)
double ccpam[4] = {4.960906e-3, 2.727436e-3, 8.392311e-4, 1.556861e-3};
double ccim = 8.978749e-2;

/*
 * calculation of the matrix of general precession from deq1 to deq2,
 * the precession angles (dzeta,dzett,dthet) are computed from the
 * constants (dc1-dc9) corresponding to the definitions in the explanatory
 * supplement to the astronomical ephemeris (1961, p.30f)
 */
static void gpre(double deq1,double deq2, double dprema[3][3])
{
    double dcsar=4.848136812e-6,dc1900=1900.0,dc1m2=0.01;
    double dc1=2304.25,dc2=1.396,dc3=0.302,dc4=0.018,dc5=0.791;
    double dc6=2004.683,dc7 = -0.853,dc8 = -0.426,dc9 = -0.042;
    double dt0,dt,dts,dtc,dzeta,dzett,dthet,dszeta,dczeta,dszett;
    double dczett,dsthet,dcthet,da,db,dc,dd;

    dt0=(deq1-dc1900)*dc1m2;
    dt=(deq2-deq1)*dc1m2;
    dts=dt*dt;
    dtc=dts*dt;
    dzeta=((dc1+dc2*dt0)*dt+dc3*dts+dc4*dtc)*dcsar;
    dzett=dzeta +dc5*dts*dcsar;
    dthet=((dc6+dc7*dt0)*dt+dc8*dts+dc9*dtc)*dcsar;
    dszeta=sin(dzeta);
    dczeta=cos(dzeta);
    dszett=sin(dzett);
    dczett=cos(dzett);
    dsthet=sin(dthet);
    dcthet=cos(dthet);
    da=dszeta*dszett;
    db=dczeta*dszett;
    dc=dszeta*dczett;
    dd=dczeta*dczett;
    dprema[0][0] = dd*dcthet-da;
    dprema[1][0] = -dc*dcthet-db;
    dprema[2][0] = -dsthet*dczett;
    dprema[0][1] = db*dcthet+dc;
    dprema[1][1] = -da*dcthet+dd;
    dprema[2][1] = -dsthet*dszett;
    dprema[0][2] = dczeta*dsthet;
    dprema[1][2] = -dszeta*dsthet;
    dprema[2][2] = dcthet;
}

/* p.stumpf, a&a sup.,41,1,1980
 *
 * calculation of the heliocentric and barycentric velocity components
 * of the earth, the largest deviations from the jpl-de96 are 42 cm/s
 * for both heliocentric and barycentric velocity components.
 *
 * given dje = julian ephemeris date
 * deq = epoch of mean equator and mean equinox of dvelh
 * and dvelb. if deq=0, both velocity vectors
 * are referred to mean equator and mean equinox of dje.
 *
 * result dvelh(k)=heliocentric, dvelb(k)=barycentric velocity
 * components, (k=1,2,3 - dx/dt,dy/dt,dz/dt. units=a.u./s)

 * the common /barxyz/ contains those intermediate results of
 * barvel which can be used to compute the heliocentric and
 * barycentric coordinates (subroutine barcor).
 */
static void barvel1(double dje, double deq, double *dvelh, double *dvelb, double *er)
{
    double sn[4];
    double dt,dtsq,dlocal,deps,a,pertl,pertld,pertr,pertrd,cosa,sina;
    double esq,dparam,param,twoe,twog,phi,f,sinf,cosf,phid,psid,drd;
    double drld,dtl,pertp,pertpd,tl,dml,b;
    double dxhd,dyhd,dzhd,dxbd,dybd,dzbd,plon,pomg,pecc;
    double dzahd,dyahd,dzabd,dyabd,deqdat;
    double e,g;
    int  k,n,k1,k3,k5,k9,k13;

    // control-parameter ideq, and time-arguments
    ideq = deq;
    dt = (dje - dct0)/dcjul;
    dtsq = dt*dt;

    // values of all elements for the instant dje
    for (k=0; k<8; k++)
    {
        dlocal=dmod(dcfel[k][0]+dt*dcfel[k][1]+dtsq*dcfel[k][2],dc2pi);
        if(k==0) dml=dlocal;
        k1=k-1;
        if(k) forbel[k1]=dlocal;
    }
    deps=dmod(dceps[0]+dt*dceps[1]+dtsq*dceps[2],dc2pi);
    for (k=0; k<17; k++)
        sorbel[k]=dmod(ccsel[k][0]+dt*ccsel[k][1]+dtsq*ccsel[k][2],dc2pi);

    //secular perturbation in longitude
    for(k=0; k<4; k++)
    {
        a=dmod(ccsec[k][1]+dt*ccsec[k][2],dc2pi);
        sn[k]=sin(a);
    }

    //periodic perturbations of the emb (earth+moon barycenter)
    pertl=ccsec[0][0]*sn[0]+ccsec[1][0]*sn[1] +
    (ccsec[2][0]+dt*ccsec3)*sn[2]+ccsec[3][0]*sn[3];
    pertld=0.0;
    pertr =0.0;
    pertrd=0.0;
    for(k=0; k<15; k++)
    {
        a=dmod(dcargs[k][0]+dt*dcargs[k][1],dc2pi);
        cosa=cos(a);
        sina=sin(a);
        pertl=pertl+ccamps[k][0]*cosa+ccamps[k][1]*sina;
        pertr=pertr+ccamps[k][2]*cosa+ccamps[k][3]*sina;
        if(k>=10) continue;
        pertld=pertld+(ccamps[k][1]*cosa-ccamps[k][0]*sina)*ccamps[k][4];
        pertrd=pertrd+(ccamps[k][3]*cosa-ccamps[k][2]*sina)*ccamps[k][4];
    }

    // elliptic part of the motion of the emb
    e = sorbel[0];
    g = forbel[0];
    esq=e*e;
    dparam=dc1-esq;
    param=dparam;
    twoe=e+e;
    twog=g+g;
    phi=twoe*((1.0-esq*0.125)*sin(g)+e*0.625*sin(twog) +
              esq*0.5416667*sin(g+twog));
    f=g+phi;
    sinf=sin(f);
    cosf=cos(f);
    dpsi=dparam/(dc1+e*cosf);
    phid=twoe*ccsgd*((1.0+esq*1.5)*cosf+e*(1.25-sinf*sinf*0.5));
    psid=ccsgd*e*sinf/sqrt(param);

    // perturbed heliocentric motion of the emb.
    d1pdro=(dc1+pertr);
    drd=d1pdro*(psid+dpsi*pertrd);
    drld=d1pdro*dpsi*(dcsld+phid+pertld);
    dtl=dmod(dml+phi+pertl,dc2pi);
    dsinls=sin(dtl);
    dcosls=cos(dtl);
    dxhd=drd*dcosls-drld*dsinls;
    dyhd=drd*dsinls+drld*dcosls;

    // influence of eccentricity, evection and variation on the
    // geocentric motion of the moon
    pertl =0.0;
    pertld=0.0;
    pertp =0.0;
    pertpd=0.0;
    for(k=0; k<3; k++)
    {
        a=dmod(dcargm[k][0]+dt*dcargm[k][1],dc2pi);
        sina=sin(a);
        cosa=cos(a);
        pertl =pertl +ccampm[k][0]*sina;
        pertld=pertld+ccampm[k][1]*cosa;
        pertp =pertp +ccampm[k][2]*cosa;
        pertpd=pertpd-ccampm[k][3]*sina;
    }

    // heliocentric motion of the earth
    tl=forbel[1]+pertl;
    sinlm=sin(tl);
    coslm=cos(tl);
    sigma=cckm/(1.0+pertp);
    a=sigma*(ccmld+pertld);
    b=sigma*pertpd;
    dxhd=dxhd+a*sinlm+b*coslm;
    dyhd=dyhd-a*coslm+b*sinlm;
    dzhd=    -sigma*ccfdi*cos(forbel[2]);

    // barycentric motion of the earth
    dxbd=dxhd*dc1mme;
    dybd=dyhd*dc1mme;
    dzbd=dzhd*dc1mme;
    for(k=0; k<4; k++)
    {
        k3=k+3;
        k1=k+1;
        k9=k+9;
        k5=k+5;
        k13=k+13;
        plon=forbel[k3];
        pomg=sorbel[k1];
        pecc=sorbel[k9];
        tl=dmod(plon+2.0*pecc*sin(plon-pomg), dc2pi);
        sinlp[k]=sin(tl);
        coslp[k]=cos(tl);
        dxbd=dxbd+ccpamv[k]*(sinlp[k]+pecc*sin(pomg));
        dybd=dybd-ccpamv[k]*(coslp[k]+pecc*cos(pomg));
        dzbd=dzbd-ccpamv[k]*sorbel[k13]*cos(plon-sorbel[k5]);
    }

    // transition to mean equator of date
    dcosep=cos(deps);
    dsinep=sin(deps);
    dyahd=dcosep*dyhd-dsinep*dzhd;
    dzahd=dsinep*dyhd+dcosep*dzhd;
    dyabd=dcosep*dybd-dsinep*dzbd;
    dzabd=dsinep*dybd+dcosep*dzbd;

    // er is the mean obliquity of the ecliptic
    *er=deps;

    if(!ideq)
    {
        dvelh[0]=dxhd;
        dvelh[1]=dyahd;
        dvelh[2]=dzahd;
        dvelb[0]=dxbd;
        dvelb[1]=dyabd;
        dvelb[2]=dzabd;
        return;
    }

    // general precession form epoch dje to deq
    deqdat=(dje-dct0-dcbes)/dctrop+1900.0;
    gpre(deqdat,deq,dprema);
    for(n=0; n<3; n++)
    {
        dvelh[n]=dxhd*dprema[0][n]+dyahd*dprema[1][n]+dzahd*dprema[2][n];
        dvelb[n]=dxbd*dprema[0][n]+dyabd*dprema[1][n]+dzabd*dprema[2][n];
    }
}

/* p.stumpff, a&a sup., 41, 1, 1980.
 *
 * calculation of heliocentric and barycentric coordinates of the earth.
 * the largest deviations from the jpl-de96 are .000011 a.u. for the
 * heliocentric, and .000046 a.u. for the barycentric coordinates.
 *
 * given the data in the common /barxyz/. they must be generated by calling
 * subroutine barvel immediately before this call of subroutine barcor.
 * result dcorh(k)=heliocentric, dcorb(k)=barycentric coordinates.
 * mean equator and mean equinox as in barvel.
 * (k=1,2,3 - x,y,z ; units= a.u.)
 */
static void barcor(double *dcorh, double *dcorb)
{
    double dr,flat,flatm,a,b,dxh,dyh,dzh,dxb,dyb,dzb,dyah,dzah,dyab,dzab;
    int  n,k,k1,k3,k5,k9,k13;

    // heliocentric coordinates of the earth - barvel
    dr=dpsi*d1pdro;
    flatm=ccim*sin(forbel[2]);
    a=sigma*cos(flatm);
    dxh=dr*dcosls - a*coslm;
    dyh=dr*dsinls - a*sinlm;
    dzh=          - sigma*sin(flatm);

    // barycentric coordinates of the earth - barvel
    dxb=dxh*dc1mme;
    dyb=dyh*dc1mme;
    dzb=dzh*dc1mme;
    for(k=0; k<4; k++)
    {
        k1=k+1;
        k3=k+3;
        k5=k+5;
        k9=k+9;
        k13=k+13;
        flat=sorbel[k13]*sin(forbel[k3]-sorbel[k5]);
        a=ccpam[k]*(1.0-sorbel[k9]*cos(forbel[k3]-sorbel[k1]));
        b=a*cos(flat);
        dxb=dxb-b*coslp[k];
        dyb=dyb-b*sinlp[k];
        dzb=dzb-a*sin(flat);
    }

    // transition to mean equator of date
    dyah=dcosep*dyh-dsinep*dzh;
    dzah=dsinep*dyh+dcosep*dzh;
    dyab=dcosep*dyb-dsinep*dzb;
    dzab=dsinep*dyb+dcosep*dzb;

    if(!ideq)
    {
        dcorh[0]=dxh;
        dcorh[1]=dyah;
        dcorh[2]=dzah;
        dcorb[0]=dxb;
        dcorb[1]=dyab;
        dcorb[2]=dzab;
        return;
    }

    // general precession from epoch dje to deq
    for(n=0; n<3; n++)
    {
        dcorh[n]=dxh*dprema[0][n]+dyah*dprema[1][n]+dzah*dprema[2][n];
        dcorb[n]=dxb*dprema[0][n]+dyab*dprema[1][n]+dzab*dprema[2][n];
    }
}

/*
 * Calculate the barycentric julian day for a specified julian day,
 * in the direction (ra,dec in radians) of coords
 */
double jdtobjd(double jd, double2 d)
{
    double  er;
    double  dvelh[3],dvelb[3],dcorh[3],dcorb[3];
    barvel1(jd, 0, dvelh, dvelb, &er);
    barcor(dcorh, dcorb);

    double adt = cos(d.y)*cos(d.x);
    double bdt = tan(er)*sin(d.y)+cos(d.y)*sin(d.x);

    // Barycentric offset in seconds
    double bc  = 499.012*(adt*dcorb[0] + bdt*dcorb[1]);

    // Heliocentric offset (unused)
    //double bh  = 499.012*(adt*dcorh[0] + bdt*dcorh[1]);

    return jd + bc/86400.0;
}

/*
 * Precess coordinates (ra,dec in radians) specified at t0 to a new epoc t1
 */
double2 precess(double2 coords, double t0, double t1)
{
    double t = -(t1-t0)/100;
    double m = (1.2812323*t + 0.0003879*t*t + 0.0000101*t*t*t)*M_PI/180;
    double n = (0.5567530*t - 0.0001185*t*t - 0.0000116*t*t*t)*M_PI/180;

    double alpham = coords.x - 1.0/2.0*(m + n*sin(coords.x)*tan(coords.y));
    double deltam = coords.y - 1.0/2.0*n*cos(alpham);

    double ra = coords.x - m - n*sin(alpham)*tan(deltam);
    double dec = coords.y - n*cos(alpham);
    return (double2) {ra, dec};
}

/*
 * Calculate the Julian day for a given UT struct tm
 */
double sumday[12]= {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
double tmtojd(struct tm *t)
{
    double days = sumday[t->tm_mon] + t->tm_mday;
    double dayfrac = t->tm_hour/24.0 + t->tm_min/1440.0 + t->tm_sec/86400.0;
    int year = t->tm_year + 1900;

    // Leap year (tm_mon = 1 for feb)
    int nleap = (t->tm_year - 1)/4;
    if (t->tm_mon > 1 && (year - 4*(year/4)) == 0)
        days += 1;

    return 2415019.5 + 365*t->tm_year + nleap + days + dayfrac;
}

/*
 * Calculate the (fractional) year for a given UT struct tm
 */
double tmtoyear(struct tm *t)
{
    double days = sumday[t->tm_mon] + t->tm_mday;
    double year = t->tm_year + 1900;

    // Leap year (tm_mon = 1 for feb)
    if (t->tm_mon > 1 && (year - 4*(year/4)) == 0)
        days += 1;

    double dayfrac = t->tm_hour/24.0 + t->tm_min/1440.0 + t->tm_sec/86400.0;
    return year + (days + dayfrac)/365.0;
}

/*
 * Calculate the Terrestrial Time offset (in seconds) for a given UTC timestamp
 * TT ~= TAI + 32.184
 */
double utcttoffset(time_t ut)
{
    size_t total = sizeof(leapseconds)/sizeof(leapseconds[0]);

    for (int i = total - 2; i >= 0; i -= 2)
        if (leapseconds[i] - 2208988800 < ut)
            return 32.184 + leapseconds[i + 1];

    error("Warning: Date before 1 Jan 1972, Leap second offset is probably incorrect");
    return 32.184;
}