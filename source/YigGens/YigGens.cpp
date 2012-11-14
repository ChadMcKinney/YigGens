/*
    YigGens are a collection of UGens by Chad McKinney created in 2011 for the network music program, Yig. Under GPL license.

    UGens:

    YigCliffordN:

        xn+1 = sin(a yn) + c cos(a xn)
        yn+1 = sin(b xn) + d cos(b yn)

        more info: http://paulbourke.net/fractals/clifford/

    Clifford 3D:
        xnew = sin(a*y) - z*cos(b*x)
        ynew = z*sin(c*x) - cos(d*y)
        znew = sin(x)

        Default parameters are: a = 2.24, b = .43, c = -.65, d = -2.43

        more info: http://www.nahee.com/spanky/www/fractint/pickover.html

    To Do:
    Then make Linear and Cubic interpolation versions.
    Next make a breakpoint version. Then make other UGens:



    Peter De Jong:
        xn+1 = sin(a yn) - cos(b xn)
        yn+1 = sin(c xn) - cos(d yn)

        more info: http://paulbourke.net/fractals/peterdejong/

    Peter De Jong Variation a la Johnny Svensson:

        xn+1 = d sin(a xn) - sin(b yn)
        yn+1 = c cos(a xn) + cos(b yn)

        more info: more info: http://paulbourke.net/fractals/peterdejong/

    Martin Attractor:

        x(n+1) = y(n) - sign(x(n))*sqrt(abs(b*x(n)-c))
        y(n+1) = a - x(n)

        http://www.nahee.com/spanky/www/fractint/martin_hop_type.html

    MandelBulb:

        r = sqrt(x*x + y*y + z*z )
        theta = atan2(sqrt(x*x + y*y) , z)
        phi = atan2(y,x)

        newx = r^n * sin(theta*n) * cos(phi*n)
        newy = r^n * sin(theta*n) * sin(phi*n)
        newz = r^n * cos(theta*n)


    more info: http://www.skytopia.com/project/fractal/2mandelbulb.html#formula & http://www.treblig.org/3dbrot/mandelbulb.html

    4D Quaternion Julia Set:

    http://paulbourke.net/fractals/quatjulia/

*/

#include "SC_PlugIn.h"

static InterfaceTable *ft;

struct YigCliffordN : public Unit
{
    float xn, yn, a, b, c, d, counter;
    int initialIterations;
};

struct YigCliffordL : public Unit
{
    float xn, xnm1, yn, ynm1, a, b, c, d, counter, frac;
    int initialIterations;
};

struct YigCliffordC : public Unit
{
    float xn, xnm1, xnm2, xnm3, yn, ynm1, ynm2, ynm3, a, b, c, d, counter, frac;
    double xc0, xc1, xc2,xc3, yc0, yc1, yc2, yc3;
    int initialIterations;
};

struct YigClifford3DN : public Unit
{
    float xn, yn, zn, a, b, c, d, counter;
    int initialIterations;
};

struct YigClifford3DL : public Unit
{
    float xn, xnm1, yn, ynm1, zn, znm1, a, b, c, d, counter, frac;
    int initialIterations;
};

struct YigClifford3DC : public Unit
{
    float xn, xnm1, xnm2, xnm3, yn, ynm1, ynm2, ynm3, zn, znm1, znm2, znm3, a, b, c, d, counter, frac;
    double xc0, xc1, xc2,xc3, yc0, yc1, yc2, yc3, zc0, zc1, zc2, zc3;
    int initialIterations;
};

struct YigMandelbulbN : public Unit
{
    float xn, yn, zn, r, dr, p, xLow, xHigh, yLow, yHigh, zLow, zHigh, fx, fy, fz, counter;
    int maxIter, size, zindex, yindex, xindex, rule;
};

extern "C"
{
    void load(InterfaceTable *inTable);
    void YigCliffordN_next(YigCliffordN *unit, int inNumSamples);
    void YigCliffordN_Ctor(YigCliffordN* unit);
    void YigCliffordL_next(YigCliffordL *unit, int inNumSamples);
    void YigCliffordL_Ctor(YigCliffordL* unit);
    void YigCliffordC_next(YigCliffordC *unit, int inNumSamples);
    void YigCliffordC_Ctor(YigCliffordC* unit);
    void YigClifford3DN_next(YigClifford3DN *unit, int inNumSamples);
    void YigClifford3DN_Ctor(YigClifford3DN* unit);
    void YigClifford3DL_next(YigClifford3DL *unit, int inNumSamples);
    void YigClifford3DL_Ctor(YigClifford3DL* unit);
    void YigClifford3DC_next(YigClifford3DC *unit, int inNumSamples);
    void YigClifford3DC_Ctor(YigClifford3DC* unit);
    void YigMandelbulbN_next(YigMandelbulbN *unit, int inNumSamples);
    void YigMandelbulbN_Ctor(YigMandelbulbN* unit);
}

////////////////////////////////////////////////////////
// calc 3rd order interpolation coefs from four points
////////////////////////////////////////////////////////
static inline void ipol3Coef(
        double xnm3, double xnm2, double xnm1, double xn,
        double &c0, double &c1, double &c2, double &c3)
{
        c0 = xnm2;
        c1 = 0.5f * (xnm1 - xnm3);
        c2 = xnm3 - (2.5f * xnm2) + xnm1 + xnm1 - 0.5f * xn;
        c3 = 0.5f * (xn - xnm3) + 1.5f * (xnm2 - xnm1);
}

// do 3rd order interpolation using coefs
static inline double ipol3(float frac, double c0, double c1, double c2, double c3){
        return ((c3 * frac + c2) * frac + c1) * frac + c0;
}

/////////////////
// YigCliffordN
////////////////

void YigCliffordN_Ctor(YigCliffordN *unit)
{
    SETCALC(YigCliffordN_next);
    unit->a = ZIN0(1);
    unit->b = ZIN0(2);
    unit->c = ZIN0(3);
    unit->d = ZIN0(4);
    unit->xn = ZIN0(5);
    unit->yn = ZIN0(6);
    unit->initialIterations = 100;
    unit->counter = 0.f;
    float a = unit->a;
    float b = unit->b;
    float c = unit->c;
    float d = unit->d;
    float xn = unit->xn;
    float yn = unit->yn;

    // Update initial conditions
    for (int i = 0; i < unit->initialIterations; i++) {
        float xnew = sin(yn*b) + c*sin(xn*b);
        float ynew = sin(xn*a) + d*sin(yn*a);
        xn = xnew;
        yn = ynew;
    }

    unit->xn = xn;
    unit->yn = yn;

    YigCliffordN_next(unit, 1);
}

void YigCliffordN_next(YigCliffordN *unit, int inNumSamples)
{
    float *leftout = ZOUT(0);
    float *rightout = ZOUT(1);
    float freq = ZIN0(0);
    float a = ZIN0(1);
    float b = ZIN0(2);
    float c = ZIN0(3);
    float d = ZIN0(4);
    float xn = unit->xn;
    float yn = unit->yn;

    float counter = unit->counter;

    float samplesPerCycle;
    if(freq < unit->mRate->mSampleRate)
    {
        samplesPerCycle = unit->mRate->mSampleRate / sc_max(freq, 0.001f);
    }
    else
    {
        samplesPerCycle = 1.f;
    }

    for(int i = 0; i < inNumSamples; i++)
    {
        if(counter >= samplesPerCycle)
        {
            counter -= samplesPerCycle;

            // compute a new point using the strange attractor equations
            float xnew = sin(yn*b) + c*sin(xn*b);
            float ynew = sin(xn*a) + d*sin(yn*a);

            // save the new point
            xn = xnew;
            yn = ynew;

        }
        counter++;
        // output the current point
        ZXP(leftout) = xn / 2;
        ZXP(rightout) = yn / 2;
    }
    unit->counter = counter;
    unit->xn = xn;
    unit->yn = yn;
}

/////////////////
// YigCliffordL
////////////////

void YigCliffordL_Ctor(YigCliffordL *unit)
{
    SETCALC(YigCliffordL_next);
    unit->a = ZIN0(1);
    unit->b = ZIN0(2);
    unit->c = ZIN0(3);
    unit->d = ZIN0(4);
    unit->xn = unit->xnm1 = ZIN0(5);
    unit->yn = unit->ynm1 = ZIN0(6);
    unit->initialIterations = 100;
    unit->counter = 0.f;
    unit->frac = 0.f;
    float a = unit->a;
    float b = unit->b;
    float c = unit->c;
    float d = unit->d;
    float xn = unit->xn;
    float yn = unit->yn;

    // Update initial conditions
    for (int i = 0; i < unit->initialIterations; i++) {
        float xnew = sin(yn*b) + c*sin(xn*b);
        float ynew = sin(xn*a) + d*sin(yn*a);
        xn = xnew;
        yn = ynew;
    }

    unit->xn = xn;
    unit->yn = yn;

    YigCliffordL_next(unit, 1);
}

void YigCliffordL_next(YigCliffordL *unit, int inNumSamples)
{
    float *leftout = ZOUT(0);
    float *rightout = ZOUT(1);
    float freq = ZIN0(0);
    float a = ZIN0(1);
    float b = ZIN0(2);
    float c = ZIN0(3);
    float d = ZIN0(4);
    float xn = unit->xn;
    float yn = unit->yn;
    float xnm1 = unit->xnm1;
    float ynm1 = unit->ynm1;

    float counter = unit->counter;
    float frac = unit->frac;

    float samplesPerCycle;
    double slope;
    if(freq < unit->mRate->mSampleRate)
    {
        samplesPerCycle = unit->mRate->mSampleRate / sc_max(freq, 0.001f);
        slope = 1.f/ samplesPerCycle;
    }
    else
    {
        samplesPerCycle = 1.f;
        slope = 1.f;
    }

    double dx = xn - xnm1;
    double dy = yn - ynm1;

    for(int i = 0; i < inNumSamples; i++)
    {
        if(counter >= samplesPerCycle)
        {
            counter -= samplesPerCycle;

            xnm1 = xn;
            ynm1 = yn;

            // compute a new point using the strange attractor equations
            float xnew = sin(yn*b) + c*sin(xn*b);
            float ynew = sin(xn*a) + d*sin(yn*a);

            // save the new point
            xn = xnew;
            yn = ynew;
            frac = 0.f;
            dx = xn - xnm1;
            dy = yn - ynm1;
        }
        counter++;
        // output the current point
        ZXP(leftout) = (xn / 2) + (dx * frac);
        ZXP(rightout) = (yn / 2) + (dy * frac);
        frac += slope;
    }
    unit->xn = xn;
    unit->xnm1 = xnm1;
    unit->yn = yn;
    unit->ynm1 = ynm1;
    unit->frac = frac;
    unit->counter = counter;
}

/////////////////
// YigCliffordC
////////////////

void YigCliffordC_Ctor(YigCliffordC *unit)
{
    SETCALC(YigCliffordC_next);
    unit->a = ZIN0(1);
    unit->b = ZIN0(2);
    unit->c = ZIN0(3);
    unit->d = ZIN0(4);
    unit->xnm1 = unit->xn = ZIN0(5);
    unit->ynm1 = unit->yn = ZIN0(6);
    unit->xnm3 = unit->xnm2 = unit->xnm1;
    unit->ynm3 = unit->ynm2 = unit->ynm1;
    unit->initialIterations = 100;
    unit->counter = 0.f;
    unit->frac = 0.f;
    float a = unit->a;
    float b = unit->b;
    float c = unit->c;
    float d = unit->d;
    float xn = unit->xn;
    float yn = unit->yn;
    unit->xc3 = unit->xc2 = unit->xc1 = unit->xc0 = 0.;
    unit->yc3 = unit->yc2 = unit->yc1 = unit->yc0 = 0.;

    // Update initial conditions
    for (int i = 0; i < unit->initialIterations; i++) {
        float xnew = sin(yn*b) + c*sin(xn*b);
        float ynew = sin(xn*a) + d*sin(yn*a);
        xn = xnew;
        yn = ynew;
    }

    unit->xn = xn;
    unit->yn = yn;

    YigCliffordC_next(unit, 1);
}

void YigCliffordC_next(YigCliffordC *unit, int inNumSamples)
{
    float *leftout = ZOUT(0);
    float *rightout = ZOUT(1);
    float freq = ZIN0(0);
    float a = ZIN0(1);
    float b = ZIN0(2);
    float c = ZIN0(3);
    float d = ZIN0(4);
    float xn = unit->xn;
    float yn = unit->yn;
    float xnm1 = unit->xnm1;
    float xnm2 = unit->xnm2;
    float xnm3 = unit->xnm3;
    float ynm1 = unit->ynm1;
    float ynm2 = unit->ynm2;
    float ynm3 = unit->ynm3;
    double xc0 = unit->xc0;
    double xc1 = unit->xc1;
    double xc2 = unit->xc2;
    double xc3 = unit->xc3;
    double yc0 = unit->yc0;
    double yc1 = unit->yc1;
    double yc2 = unit->yc2;
    double yc3 = unit->yc3;

    float counter = unit->counter;
    float frac = unit->frac;

    float samplesPerCycle;
    double slope;
    if(freq < unit->mRate->mSampleRate)
    {
        samplesPerCycle = unit->mRate->mSampleRate / sc_max(freq, 0.001f);
        slope = 1.f/ samplesPerCycle;
    }
    else
    {
        samplesPerCycle = 1.f;
        slope = 1.f;
    }

    for(int i = 0; i < inNumSamples; i++)
    {
        if(counter >= samplesPerCycle)
        {
            counter -= samplesPerCycle;
            xnm3 = xnm2;
            xnm2 = xnm1;
            xnm1 = xn;
            ynm3 = ynm2;
            ynm2 = ynm1;
            ynm1 = yn;

            // compute a new point using the strange attractor equations
            float xnew = sin(yn*b) + c*sin(xn*b);
            float ynew = sin(xn*a) + d*sin(yn*a);

            // save the new point
            xn = xnew;
            yn = ynew;
            ipol3Coef(xnm3, xnm2, xnm1, xn, xc0, xc1, xc2, xc3);
            ipol3Coef(ynm3, ynm2, ynm1, yn, yc0, yc1, yc2, yc3);
            frac = 0.f;
        }
        counter++;
        // output the current point
        //ZXP(leftout) = (xn / 2) + (dx * frac);
        //ZXP(rightout) = (yn / 2) + (dy * frac);
        ZXP(leftout) = ipol3(frac, xc0, xc1, xc2, xc3);
        ZXP(rightout) = ipol3(frac, yc0, yc1, yc2, yc3);
        frac += slope;
    }
    unit->xn = xn;
    unit->xnm1 = xnm1;
    unit->xnm2 = xnm2;
    unit->xnm3 = xnm3;
    unit->yn = yn;
    unit->ynm1 = ynm1;
    unit->ynm2 = ynm2;
    unit->ynm3 = ynm3;
    unit->frac = frac;
    unit->counter = counter;
    unit->xc0 = xc0; unit->xc1 = xc1; unit->xc2 = xc2; unit->xc3 = xc3;
    unit->yc0 = yc0; unit->yc1 = yc1; unit->yc2 = yc2; unit->yc3 = yc3;
}

//////////////////
// YigClifford3DN
/////////////////

void YigClifford3DN_Ctor(YigClifford3DN *unit)
{
    SETCALC(YigClifford3DN_next);
    unit->a = ZIN0(1);
    unit->b = ZIN0(2);
    unit->c = ZIN0(3);
    unit->d = ZIN0(4);
    unit->xn = ZIN0(5);
    unit->yn = ZIN0(6);
    unit->zn = ZIN0(7);
    unit->initialIterations = 100;
    unit->counter = 0.f;
    float a = unit->a;
    float b = unit->b;
    float c = unit->c;
    float d = unit->d;
    float xn = unit->xn;
    float yn = unit->yn;
    float zn = unit->zn;

    // Update initial conditions
    for (int i = 0; i < unit->initialIterations; i++) {

        float xnew = sin(a*yn) - zn*cos(b*xn);
        float ynew = zn*sin(c*xn) - cos(d*yn);
        float znew = sin(xn);

        xn = xnew;
        yn = ynew;
        zn = znew;
    }

    unit->xn = xn;
    unit->yn = yn;
    unit->zn = zn;

    YigClifford3DN_next(unit, 1);
}

void YigClifford3DN_next(YigClifford3DN *unit, int inNumSamples)
{
    float *xout = ZOUT(0);
    float *yout = ZOUT(1);
    float *zout = ZOUT(2);
    float freq = ZIN0(0);
    float a = ZIN0(1);
    float b = ZIN0(2);
    float c = ZIN0(3);
    float d = ZIN0(4);
    float xn = unit->xn;
    float yn = unit->yn;
    float zn = unit->zn;

    float counter = unit->counter;

    float samplesPerCycle;
    if(freq < unit->mRate->mSampleRate)
    {
        samplesPerCycle = unit->mRate->mSampleRate / sc_max(freq, 0.001f);
    }
    else
    {
        samplesPerCycle = 1.f;
    }

    for(int i = 0; i < inNumSamples; i++)
    {
        if(counter >= samplesPerCycle)
        {
            counter -= samplesPerCycle;

            float xnew = sin(a*yn) - zn*cos(b*xn);
            float ynew = zn*sin(c*xn) - cos(d*yn);
            float znew = sin(xn);

            xn = xnew;
            yn = ynew;
            zn = znew;

        }
        counter++;
        // output the current point
        ZXP(xout) = xn / 2;
        ZXP(yout) = yn / 2;
        ZXP(zout) = zn / 2;
    }
    unit->counter = counter;
    unit->xn = xn;
    unit->yn = yn;
    unit->zn = zn;
}

//////////////////
// YigClifford3DL
/////////////////

void YigClifford3DL_Ctor(YigClifford3DL *unit)
{
    SETCALC(YigClifford3DL_next);
    unit->a = ZIN0(1);
    unit->b = ZIN0(2);
    unit->c = ZIN0(3);
    unit->d = ZIN0(4);
    unit->xn = unit->xnm1 = ZIN0(5);
    unit->yn = unit->ynm1 = ZIN0(6);
    unit->zn = unit->znm1 = ZIN0(7);
    unit->initialIterations = 100;
    unit->counter = 0.f;
    unit->frac = 0.f;
    float a = unit->a;
    float b = unit->b;
    float c = unit->c;
    float d = unit->d;
    float xn = unit->xn;
    float yn = unit->yn;
    float zn = unit->zn;

    // Update initial conditions
    for (int i = 0; i < unit->initialIterations; i++) {

        float xnew = sin(a*yn) - zn*cos(b*xn);
        float ynew = zn*sin(c*xn) - cos(d*yn);
        float znew = sin(xn);

        xn = xnew;
        yn = ynew;
        zn = znew;
    }

    unit->xn = xn;
    unit->yn = yn;
    unit->zn = zn;

    YigClifford3DL_next(unit, 1);
}

void YigClifford3DL_next(YigClifford3DL *unit, int inNumSamples)
{
    float *xout = ZOUT(0);
    float *yout = ZOUT(1);
    float *zout = ZOUT(2);
    float freq = ZIN0(0);
    float a = ZIN0(1);
    float b = ZIN0(2);
    float c = ZIN0(3);
    float d = ZIN0(4);
    float xn = unit->xn;
    float yn = unit->yn;
    float zn = unit->zn;
    float xnm1 = unit->xnm1;
    float ynm1 = unit->ynm1;
    float znm1 = unit->znm1;

    float counter = unit->counter;
    float frac = unit->frac;

    float samplesPerCycle;
    double slope;
    if(freq < unit->mRate->mSampleRate)
    {
        samplesPerCycle = unit->mRate->mSampleRate / sc_max(freq, 0.001f);
        slope = 1.f/ samplesPerCycle;
    }
    else
    {
        samplesPerCycle = 1.f;
        slope = 1.f;
    }

    double dx = xn - xnm1;
    double dy = yn - ynm1;
    double dz = zn - znm1;

    for(int i = 0; i < inNumSamples; i++)
    {
        if(counter >= samplesPerCycle)
        {
            counter -= samplesPerCycle;

            xnm1 = xn;
            ynm1 = yn;
            znm1 = zn;

            // compute a new point using the strange attractor equations
            float xnew = sin(a*yn) - zn*cos(b*xn);
            float ynew = zn*sin(c*xn) - cos(d*yn);
            float znew = sin(xn);

            xn = xnew;
            yn = ynew;
            zn = znew;

            frac = 0.f;
            dx = xn - xnm1;
            dy = yn - ynm1;
            dz = zn - znm1;
        }
        counter++;
        // output the current point
        ZXP(xout) = (xn / 2) + (dx * frac);
        ZXP(yout) = (yn / 2) + (dy * frac);
        ZXP(zout) = (zn / 2) + (dz * frac);
        frac += slope;
    }
    unit->xn = xn;
    unit->xnm1 = xnm1;
    unit->yn = yn;
    unit->ynm1 = ynm1;
    unit->zn = zn;
    unit->znm1 = znm1;
    unit->frac = frac;
    unit->counter = counter;
}

//////////////////
// YigClifford3DC
/////////////////

void YigClifford3DC_Ctor(YigClifford3DC *unit)
{
    SETCALC(YigClifford3DC_next);
    unit->a = ZIN0(1);
    unit->b = ZIN0(2);
    unit->c = ZIN0(3);
    unit->d = ZIN0(4);
    unit->xnm1 = unit->xn = ZIN0(5);
    unit->ynm1 = unit->yn = ZIN0(6);
    unit->znm1 = unit->zn = ZIN0(7);
    unit->xnm3 = unit->xnm2 = unit->xnm1;
    unit->ynm3 = unit->ynm2 = unit->ynm1;
    unit->znm3 = unit->znm2 = unit->znm1;
    unit->initialIterations = 100;
    unit->counter = 0.f;
    unit->frac = 0.f;
    float a = unit->a;
    float b = unit->b;
    float c = unit->c;
    float d = unit->d;
    float xn = unit->xn;
    float yn = unit->yn;
    float zn = unit->zn;
    unit->xc3 = unit->xc2 = unit->xc1 = unit->xc0 = 0.;
    unit->yc3 = unit->yc2 = unit->yc1 = unit->yc0 = 0.;
    unit->zc3 = unit->zc2 = unit->zc1 = unit->zc0 = 0.;

    // Update initial conditions
    for (int i = 0; i < unit->initialIterations; i++) {

        float xnew = sin(a*yn) - zn*cos(b*xn);
        float ynew = zn*sin(c*xn) - cos(d*yn);
        float znew = sin(xn);

        xn = xnew;
        yn = ynew;
        zn = znew;
    }

    unit->xn = xn;
    unit->yn = yn;
    unit->zn = zn;

    YigClifford3DC_next(unit, 1);
}

void YigClifford3DC_next(YigClifford3DC *unit, int inNumSamples)
{
    float *xout = ZOUT(0);
    float *yout = ZOUT(1);
    float *zout = ZOUT(2);
    float freq = ZIN0(0);
    float a = ZIN0(1);
    float b = ZIN0(2);
    float c = ZIN0(3);
    float d = ZIN0(4);
    float xn = unit->xn;
    float yn = unit->yn;
    float zn = unit->zn;
    float xnm1 = unit->xnm1;
    float xnm2 = unit->xnm2;
    float xnm3 = unit->xnm3;
    float ynm1 = unit->ynm1;
    float ynm2 = unit->ynm2;
    float ynm3 = unit->ynm3;
    float znm1 = unit->znm1;
    float znm2 = unit->znm2;
    float znm3 = unit->znm3;
    double xc0 = unit->xc0;
    double xc1 = unit->xc1;
    double xc2 = unit->xc2;
    double xc3 = unit->xc3;
    double yc0 = unit->yc0;
    double yc1 = unit->yc1;
    double yc2 = unit->yc2;
    double yc3 = unit->yc3;
    double zc0 = unit->zc0;
    double zc1 = unit->zc1;
    double zc2 = unit->zc2;
    double zc3 = unit->zc3;

    float counter = unit->counter;
    float frac = unit->frac;

    float samplesPerCycle;
    double slope;
    if(freq < unit->mRate->mSampleRate)
    {
        samplesPerCycle = unit->mRate->mSampleRate / sc_max(freq, 0.001f);
        slope = 1.f/ samplesPerCycle;
    }
    else
    {
        samplesPerCycle = 1.f;
        slope = 1.f;
    }

    for(int i = 0; i < inNumSamples; i++)
    {
        if(counter >= samplesPerCycle)
        {
            counter -= samplesPerCycle;
            xnm3 = xnm2;
            xnm2 = xnm1;
            xnm1 = xn;
            ynm3 = ynm2;
            ynm2 = ynm1;
            ynm1 = yn;
            znm3 = znm2;
            znm2 = znm1;
            znm1 = zn;

            // compute a new point using the strange attractor equations
            float xnew = sin(a*yn) - zn*cos(b*xn);
            float ynew = zn*sin(c*xn) - cos(d*yn);
            float znew = sin(xn);

            xn = xnew;
            yn = ynew;
            zn = znew;

            ipol3Coef(xnm3, xnm2, xnm1, xn, xc0, xc1, xc2, xc3);
            ipol3Coef(ynm3, ynm2, ynm1, yn, yc0, yc1, yc2, yc3);
            ipol3Coef(znm3, znm2, znm1, zn, zc0, zc1, zc2, zc3);
            frac = 0.f;
        }
        counter++;
        // output the current point
        //ZXP(leftout) = (xn / 2) + (dx * frac);
        //ZXP(rightout) = (yn / 2) + (dy * frac);
        ZXP(xout) = ipol3(frac, xc0, xc1, xc2, xc3);
        ZXP(yout) = ipol3(frac, yc0, yc1, yc2, yc3);
        ZXP(zout) = ipol3(frac, zc0, zc1, zc2, zc3);
        frac += slope;
    }
    unit->xn = xn;
    unit->xnm1 = xnm1;
    unit->xnm2 = xnm2;
    unit->xnm3 = xnm3;
    unit->yn = yn;
    unit->ynm1 = ynm1;
    unit->ynm2 = ynm2;
    unit->ynm3 = ynm3;
    unit->zn = zn;
    unit->znm1 = znm1;
    unit->znm2 = znm2;
    unit->znm3 = znm3;
    unit->frac = frac;
    unit->counter = counter;
    unit->xc0 = xc0; unit->xc1 = xc1; unit->xc2 = xc2; unit->xc3 = xc3;
    unit->yc0 = yc0; unit->yc1 = yc1; unit->yc2 = yc2; unit->yc3 = yc3;
    unit->zc0 = zc0; unit->zc1 = zc1; unit->zc2 = zc2; unit->zc3 = zc3;
}

//////////////////
// YigMandelbulbN
/////////////////

float checkRange(float low, float high, int size, int off)
{
    return low+((high - low) / (float) size) * (float) off;
}

void YigMandelbulbN_Ctor(YigMandelbulbN *unit)
{
    SETCALC(YigMandelbulbN_next);
    unit->p = ZIN0(1);
    unit->size = (int) ZIN0(2);
    unit->rule = (int) ZIN0(3);
    unit->maxIter = (int) ZIN0(4);
    unit->xn = 0;
    unit->yn = 0;
    unit->zn = 0.9;
    unit->r = 0.f;
    unit->dr = 1.0f;
    unit->counter = 0.f;
    unit->xLow = unit->yLow = unit->zLow = -1.2f;
    unit->xHigh = unit->yHigh = unit->zHigh = 1.2f;
    unit->xindex = unit->size/2;
    unit->yindex = unit->size/2;
    unit->zindex = unit->size/2;
    unit->fx = checkRange(unit->xLow, unit->xHigh, unit->size, unit->xindex);
    unit->fy = checkRange(unit->yLow, unit->yHigh, unit->size, unit->yindex);
    unit->fz = checkRange(unit->zLow, unit->zHigh, unit->size, unit->zindex);

    YigMandelbulbN_next(unit, 1);
}

void YigMandelbulbN_next(YigMandelbulbN *unit, int inNumSamples)
{
    float *xout = ZOUT(0);
    float *yout = ZOUT(1);
    float *zout = ZOUT(2);
    float freq = ZIN0(0);
    float p = ZIN0(1);
    int size = (int) ZIN0(2);
    int rule = (int) ZIN0(3);
    int maxIter = (int) ZIN0(4);
    float xn = unit->xn;
    float yn = unit->yn;
    float zn = unit->zn;
    int xindex = unit->xindex;
    int yindex = unit->yindex;
    int zindex = unit->zindex;
    float xLow = unit->xLow;
    float yLow = unit->yLow;
    float zLow = unit->zLow;
    float xHigh = unit->xHigh;
    float yHigh = unit->yHigh;
    float zHigh= unit->zHigh;
    float fx = unit->fx;
    float fy = unit->fy;
    float fz = unit->fz;

    float counter = unit->counter;

    if(rule != unit->rule)
    {
        xindex = size - xindex;
        yindex = size - yindex;
        zindex = size - zindex;
        unit->rule = rule;
    }

    float samplesPerCycle;
    if(freq < unit->mRate->mSampleRate)
    {
        samplesPerCycle = unit->mRate->mSampleRate / sc_max(freq, 0.001f);
    }
    else
    {
        samplesPerCycle = 1.f;
    }

    for(int iter = 0; iter < inNumSamples; iter++)
    {
        if(counter >= samplesPerCycle)
        {

            counter -= samplesPerCycle;

            // iterate in 3 dimensions
            xindex++;
            fx = checkRange(xLow, xHigh, size, xindex);

            float x, y, z, xnew, ynew, znew, r, theta, phi;
            int pointVal;
            for(pointVal = 0., x = xn, y = yn, z = zn; (pointVal < maxIter) && (x * x + y * y + z * z) < 2.0; pointVal++)
            {
                /*
                float rpow;

                r = sqrt(x * x + y * y + z * z );
                theta = atan2(sqrt(x * x + y * y) , z);
                phi = atan2(y, x);

                rpow = pow(r, p);
                xnew = rpow * sin(theta * p) * cos(phi * p);
                ynew = rpow * sin(theta * p) * sin(phi * p);
                znew = rpow * cos(theta * p);

                x = xnew + fx;
                y = ynew + fy;
                z = znew + fz;*/
                r = sqrt(x * x + y * y + z * z );
                theta = atan2(y, x);
                phi = sin(z / r);

                r = pow(r, p);
                theta *= p;
                phi *= p;
                int mul;
                float angle;

                // Iterate through all possible rule combinations
                for(int j = 15; j > 0; j--)
                {
                    bool ruleToggle = (rule >> (16 - j)) & 1;
                    //qDebug() << ruleToggle;
                    switch(j)
                    {
                    case 0:

                        break;

                    case 1:
                        if(ruleToggle)
                        {
                            xnew = mul * sin(angle);
                        }
                        else
                        {
                            xnew = mul * cos(angle);
                        }
                        break;

                    case 2:
                        if(ruleToggle)
                        {
                            angle = theta;
                        }
                        else
                        {
                            angle = phi;
                        }
                        break;

                    case 3:
                        if(ruleToggle)
                        {
                            mul = 1;
                        }
                        else
                        {
                            mul = -1;
                        }
                        break;

                    case 4:
                        if(ruleToggle)
                        {
                            ynew *= mul * sin(angle);
                        }
                        else
                        {
                            ynew *= mul * cos(angle);
                        }
                        break;

                    case 5:
                        if(ruleToggle)
                        {
                            angle = theta;
                        }
                        else
                        {
                            angle = phi;
                        }
                        break;

                    case 6:
                        if(ruleToggle)
                        {
                            mul = 1;
                        }
                        else
                        {
                            mul = -1;
                        }
                        break;

                    case 7:
                        if(ruleToggle)
                        {
                            ynew = mul * sin(angle);
                        }
                        else
                        {
                            ynew = mul * cos(angle);
                        }
                        break;

                    case 8:
                        if(ruleToggle)
                        {
                            angle = theta;
                        }
                        else
                        {
                            angle = phi;
                        }
                        break;

                    case 9:
                        if(ruleToggle)
                        {
                            mul = -1;
                        }
                        else
                        {
                            mul = 1;
                        }
                        break;

                    case 10:
                        if(ruleToggle)
                        {
                            znew *= mul * sin(angle);
                        }
                        else
                        {
                            znew *= mul * cos(angle);
                        }
                        break;

                    case 11:
                        if(ruleToggle)
                        {
                            angle = theta;
                        }
                        else
                        {
                            angle = phi;
                        }
                        break;

                    case 12:
                        if(ruleToggle)
                        {
                            mul = -1;
                        }
                        else
                        {
                            mul = 1;
                        }
                        break;

                    case 13:
                        if(ruleToggle)
                        {
                            znew = mul * sin(angle);
                        }
                        else
                        {
                            znew = mul * cos(angle);
                        }
                        break;

                    case 14:
                        if(ruleToggle)
                        {
                            angle = theta;
                        }
                        else
                        {
                            angle = phi;
                        }
                        break;

                    case 15:
                        if(ruleToggle)
                        {
                            mul = -1;
                        }
                        else
                        {
                            mul = 1;
                        }
                        break;


                    }
                }
                x = (xnew * r) + fx;
                y = (ynew * r) + fy;
                z = (znew * r) + fz;

            }

            if( pointVal >= (maxIter - 1))
            {
                xn = x;
                yn = y;
                zn = z;
            }

            if(xindex >= size)
            {
                xindex = 0;
                yindex++;
                fy = checkRange(yLow, yHigh, size, yindex);
                //printf("yindex++");
            }

            if(yindex >= size)
            {
                yindex = 0;
                zindex++;
                fz = checkRange(zLow, zHigh, size, zindex);
                //printf("zindex++");
            }

            if(zindex >= size)
            {
                zindex = 0;
                //fx = 0;
                //fy = 0;
                //fz = 0;
                //printf("restart");
            }

        }
        counter++;
        // constrain and output the current point
        ZXP(xout) = xn;
        ZXP(yout) = yn;
        ZXP(zout) = zn;
    }
    unit->counter = counter;
    unit->xn = xn;
    unit->yn = yn;
    unit->zn = zn;
    unit->size = size;
    unit->xindex = xindex;
    unit->yindex = yindex;
    unit->zindex = zindex;
    unit->fx = fx;
    unit->fy = fy;
    unit->fz = fz;
    unit->rule = rule;
}

PluginLoad(Yig)
{
    ft = inTable;
    DefineSimpleUnit(YigCliffordN);
    DefineSimpleUnit(YigCliffordL);
    DefineSimpleUnit(YigCliffordC);
    DefineSimpleUnit(YigClifford3DN);
    DefineSimpleUnit(YigClifford3DL);
    DefineSimpleUnit(YigClifford3DC);
    DefineSimpleUnit(YigMandelbulbN);
}

