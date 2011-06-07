/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <rtengine.h>
#include <improcfun.h>
#include <curves.h>
#include <math.h>
#include <colorclip.h>
#include <gauss.h>
#include <bilateral2.h>
#include <minmax.h>
#include <mytime.h>
#include <glibmm.h>
#include <iccstore.h>
#include <impulse_denoise.h>
#include <utils.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace rtengine {

using namespace procparams;

#undef MAXVAL
#undef CMAXVAL
#undef MAXL
#undef MAX
#undef MIN
#undef ABS
#undef CLIP
#undef CLIPS
#undef CLIPC
#undef CLIPTO

#define MAXVAL  0xffff
#define CMAXVAL 0xffff
#define MAXL 	0xffff
#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))
#define ABS(a) ((a)<0?-(a):(a))
#define CLIP(a) ((a)>0?((a)<CMAXVAL?(a):CMAXVAL):0)
#define CLIPS(a) ((a)>-32768?((a)<32767?(a):32767):-32768)
#define CLIPC(a) ((a)>-32000?((a)<32000?(a):32000):-32000)
#define CLIPTO(a,b,c) ((a)>(b)?((a)<(c)?(a):(c)):(b))

extern const Settings* settings;

int* ImProcFunctions::cacheL = 0;
int* ImProcFunctions::cachea = 0;
int* ImProcFunctions::cacheb = 0;
int* ImProcFunctions::xcache = 0;
int* ImProcFunctions::ycache = 0;
int* ImProcFunctions::zcache = 0;
unsigned short* ImProcFunctions::gamma2curve = 0;

void ImProcFunctions::initCache () {

    const int maxindex = 2*65536;
    cacheL = new int[maxindex];
    cachea = new int[maxindex];
    cacheb = new int[maxindex];
    gamma2curve = new unsigned short[65536];

    int threshold = (int)(0.008856*CMAXVAL);
    for (int i=0; i<maxindex; i++)
        if (i>threshold) {
            cacheL[i] = (int)round(655.35 * (116.0 * exp(1.0/3.0 * log((double)i / CMAXVAL)) - 16.0));
            cachea[i] = (int)round(32768.0 * 500.0 * exp(1.0/3.0 * log((double)i / CMAXVAL)));
            cacheb[i] = (int)round(32768.0 * 200.0 * exp(1.0/3.0 * log((double)i / CMAXVAL)));
        }
        else {
            cacheL[i] = (int)round(9033.0 * (double)i / 1000.0); // assuming CMAXVAL = 65535
            cachea[i] = (int)round(32768.0 * 500.0 * (7.787*i/CMAXVAL+16.0/116.0));
            cacheb[i] = (int)round(32768.0 * 200.0 * (7.787*i/CMAXVAL+16.0/116.0));
        }

    double fY;
    ycache = new int[0x10000];
    for (int i=0; i<0x10000; i++)
        ycache[i] = (int)round(65536.0 * ((fY=((double)i/655.35+16)/116) > 2.0689655172413793e-1 ? fY*fY*fY : 1.107056459879453852e-3*(double)i/655.35));
    for (int i=0; i<0x10000; i++)
        ycache[i] = CLIP(ycache[i]);
    xcache = new int[369621];
    for (int i=-141556; i<228064; i++)
        xcache[i+141556] = (int)round(65536.0 * (i > 15728 ? ((double)i/76021)*((double)i/76021)*((double)i/76021)*0.96422 : (1.2841854934601665e-1*(double)i/76021-1.7712903358071262e-2)*0.96422));
    for (int i=0; i<369620; i++)
        xcache[i] = CLIP(xcache[i]);
    zcache = new int[825747];
    for (int i=-369619; i<456128; i++)
        zcache[i+369619] = (int)round(65536.0 * (i > 15728 ? ((double)i/76021)*((double)i/76021)*((double)i/76021)*0.82521 : (1.2841854934601665e-1*(double)i/76021-1.7712903358071262e-2)*0.82521));
    for (int i=0; i<825747; i++)
        zcache[i] = CLIP(zcache[i]);

	for (int i=0; i<65536; i++) {
		int g = (int)(CurveFactory::gamma2(i/65535.0) * 65535.0);
		gamma2curve[i] = CLIP(g);
	}
}

void ImProcFunctions::cleanupCache () {

	delete [] cacheL;
	delete [] cachea;
	delete [] cacheb;
	delete [] xcache;
	delete [] ycache;
	delete [] zcache;
	delete [] gamma2curve;
}

ImProcFunctions::~ImProcFunctions () {

	if (monitorTransform!=NULL)
		cmsDeleteTransform (monitorTransform);
}

void ImProcFunctions::setScale (double iscale) {
	scale = iscale;
}

void ImProcFunctions::firstAnalysis_ (Image16* original, const TMatrix &wprof, unsigned int* histogram, int* chroma_radius, int row_from, int row_to) {

    int toxyz[3][3];
    toxyz[0][0] = round(32768.0 * wprof[0][0] / 0.96422); 
    toxyz[1][0] = round(32768.0 * wprof[1][0] / 0.96422); 
    toxyz[2][0] = round(32768.0 * wprof[2][0] / 0.96422); 
    toxyz[0][1] = round(32768.0 * wprof[0][1]); 
    toxyz[1][1] = round(32768.0 * wprof[1][1]); 
    toxyz[2][1] = round(32768.0 * wprof[2][1]); 
    toxyz[0][2] = round(32768.0 * wprof[0][2] / 0.82521); 
    toxyz[1][2] = round(32768.0 * wprof[1][2] / 0.82521); 
    toxyz[2][2] = round(32768.0 * wprof[2][2] / 0.82521); 

	lumimul[0] = wprof[0][1];
	lumimul[1] = wprof[1][1];
	lumimul[2] = wprof[2][1];
	
    int W = original->width;
    int cradius = 1;
    for (int i=row_from; i<row_to; i++) {
        for (int j=0; j<W; j++) {
      
            int r = original->r[i][j];
            int g = original->g[i][j];
            int b = original->b[i][j];

            int x = (toxyz[0][0] * r + toxyz[1][0] * g + toxyz[2][0] * b) >> 15;
            int y = (toxyz[0][1] * r + toxyz[1][1] * g + toxyz[2][1] * b) >> 15;
            int z = (toxyz[0][2] * r + toxyz[1][2] * g + toxyz[2][2] * b) >> 15;

            x = CLIPTO(x,0,2*65536-1);
            y = CLIPTO(y,0,2*65536-1);
            z = CLIPTO(z,0,2*65536-1);

            int oa = cachea[x] - cachea[y];
            int ob = cacheb[y] - cacheb[z];

            if (oa<0) oa = -oa;
            if (ob<0) ob = -ob;

            if (oa > cradius)
                cradius = oa;
            if (ob > cradius)
                cradius = ob;

            if (histogram) {
                int hval = CLIP(y); //(306 * original->r[i][j] + 601 * original->g[i][j] + 117 * original->b[i][j]) >> 10;
                histogram[hval]++;
            }
        }
    }
    *chroma_radius = cradius;
}

void ImProcFunctions::firstAnalysis (Image16* original, const ProcParams* params, unsigned int* histogram, double gamma) {

	// set up monitor transform
	TMatrix wprof = iccStore->workingSpaceMatrix (params->icm.working);
	if (monitorTransform)
		cmsDeleteTransform (monitorTransform);
	monitorTransform = NULL;
	cmsHPROFILE monitor = iccStore->getProfile ("file:"+settings->monitorProfile);
	if (monitor) {
        cmsHPROFILE iprof = iccStore->getXYZProfile ();       
		cmsHPROFILE oprof = iccStore->getProfile (params->icm.output);
		if (!oprof)
			oprof = iccStore->getsRGBProfile ();
        lcmsMutex->lock ();
		monitorTransform = cmsCreateTransform (iprof, TYPE_RGB_16, monitor, TYPE_RGB_8, settings->colorimetricIntent, 0);
        lcmsMutex->unlock ();
	}

	// calculate chroma radius and histogram of the y channel needed for exposure curve calculation

#ifdef _OPENMP
    int T = omp_get_max_threads();
#else
    int T = 1;
#endif

	int* cr = new int [T];
    unsigned int** hist = new unsigned int* [T];
    for (int i=0; i<T; i++) {
		cr[i] = 0;
		hist[i] = new unsigned int[65536];
		memset (hist[i], 0, 65536*sizeof(int));
    }

    int H = original->height;
#ifdef _OPENMP
	#pragma omp parallel if (multiThread)
    {
		int tid = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		int blk = H/nthreads;

		if (tid<nthreads-1)
			firstAnalysis_ (original, wprof, hist[tid], &cr[tid], tid*blk, (tid+1)*blk);
		else
			firstAnalysis_ (original, wprof, hist[tid], &cr[tid], tid*blk, H);
    }
#else
    firstAnalysis_ (original, wprof, hist[0], &cr[0], 0, original->height);
#endif
    chroma_radius = cr[0];
    for (int i=0; i<T; i++)
    	if (cr[i]>chroma_radius)
    		chroma_radius = cr[i];
 
	memset (histogram, 0, 65536*sizeof(int));
    for (int i=0; i<65536; i++)
    	for (int j=0; j<T; j++)
    		histogram[i] += hist[j][i];

    chroma_scale = 100;    // 32768*32768 / (3*chroma_radius);
	//printf ("chroma_radius= %d   chroma_scale= %d\n",chroma_radius,chroma_scale);

    delete [] cr;
    for (int i=0; i<T; i++)
    	delete [] hist[i];
    delete [] hist;
}

void ImProcFunctions::rgbProc (Image16* working, LabImage* lab, float* hltonecurve, float* shtonecurve, int* tonecurve, SHMap* shmap, /*float defmul,*/ int sat) {

    int h_th, s_th;
    if (shmap) {
        h_th = shmap->max - params->sh.htonalwidth * (shmap->max - shmap->avg) / 100;
        s_th = params->sh.stonalwidth * (shmap->avg - shmap->min) / 100;
    }

    bool processSH  = params->sh.enabled && shmap!=NULL && (params->sh.highlights>0 || params->sh.shadows>0);
    bool processLCE = params->sh.enabled && shmap!=NULL && params->sh.localcontrast>0;
    double lceamount = params->sh.localcontrast / 200.0;

    TMatrix wprof = iccStore->workingSpaceMatrix (params->icm.working);
    int toxyz[3][3] = {
        {
        	floor(32768.0 * wprof[0][0] / 0.96422),
        	floor(32768.0 * wprof[0][1]),
        	floor(32768.0 * wprof[0][2] / 0.82521)
        },{
			floor(32768.0 * wprof[1][0] / 0.96422),
			floor(32768.0 * wprof[1][1]),
			floor(32768.0 * wprof[1][2] / 0.82521)
        },{
			floor(32768.0 * wprof[2][0] / 0.96422),
			floor(32768.0 * wprof[2][1]),
			floor(32768.0 * wprof[2][2] / 0.82521)
        }
    };

    bool mixchannels = params->chmixer.red[0]!=100 || params->chmixer.red[1]!=0 || params->chmixer.red[2]!=0 || params->chmixer.green[0]!=0 || params->chmixer.green[1]!=100 || params->chmixer.green[2]!=0 || params->chmixer.blue[0]!=0 || params->chmixer.blue[1]!=0 || params->chmixer.blue[2]!=100;

    int mapval;
    double factor;
    int tW = working->width;
    int tH = working->height;
    int r, g, b;
	float h, s, v;
	double pi = M_PI;
	FlatCurve* hCurve;
	FlatCurve* sCurve;
	FlatCurve* vCurve;
	
	float* cossq = new float [8093];
	for (int i=0; i<8093; i++) 
		cossq[i] = SQR(cos(pi*(float)i/16384));
	
	FlatCurveType hCurveType = (FlatCurveType)params->hsvequalizer.hcurve.at(0);
	FlatCurveType sCurveType = (FlatCurveType)params->hsvequalizer.scurve.at(0);
	FlatCurveType vCurveType = (FlatCurveType)params->hsvequalizer.vcurve.at(0);
	bool hCurveEnabled = hCurveType > FCT_Linear;
	bool sCurveEnabled = sCurveType > FCT_Linear;
	bool vCurveEnabled = vCurveType > FCT_Linear;

	// TODO: We should create a 'skip' value like for CurveFactory::complexsgnCurve (rtengine/curves.cc)
	if (hCurveEnabled) hCurve = new FlatCurve(params->hsvequalizer.hcurve);
	if (sCurveEnabled) sCurve = new FlatCurve(params->hsvequalizer.scurve);
	if (vCurveEnabled) vCurve = new FlatCurve(params->hsvequalizer.vcurve);

#pragma omp parallel for  private(r, g, b,factor,mapval,h,s,v) if (multiThread)
    for (int i=0; i<tH; i++) {

        for (int j=0; j<tW; j++) {

            r = working->r[i][j];
            g = working->g[i][j];
            b = working->b[i][j];

            if (mixchannels) {
                int newr = (r*params->chmixer.red[0]   + g*params->chmixer.red[1]   + b*params->chmixer.red[2]) / 100;
                int newg = (r*params->chmixer.green[0] + g*params->chmixer.green[1] + b*params->chmixer.green[2]) / 100;
                int newb = (r*params->chmixer.blue[0]  + g*params->chmixer.blue[1]  + b*params->chmixer.blue[2]) / 100;
                r = CLIP(newr);
                g = CLIP(newg);
                b = CLIP(newb);
            }

            if (processSH || processLCE) {
                mapval = shmap->map[i][j];
                factor = 1.0;
                
                if (processSH) {
                    if (mapval > h_th) 
                        factor = (h_th + (100.0 - params->sh.highlights) * (mapval - h_th) / 100.0) / mapval; 
                    else if (mapval < s_th) 
                        factor = (s_th - (100.0 - params->sh.shadows) * (s_th - mapval) / 100.0) / mapval; 
                }
                if (processLCE) {
                    double sub = lceamount*(mapval-factor*(r*lumimul[0] + g*lumimul[1] + b*lumimul[2]));
                    r = CLIP((int)(factor*r-sub));
                    g = CLIP((int)(factor*g-sub));
                    b = CLIP((int)(factor*b-sub));
                }
                else {
                    r = CLIP((int)(factor*r));
                    g = CLIP((int)(factor*g));
                    b = CLIP((int)(factor*b));
                }
            }


			float tonefactor=(hltonecurve[r]+hltonecurve[g]+hltonecurve[b])/3;

			r = (r*tonefactor);
			g = (g*tonefactor);
			b = (b*tonefactor);
			
			//shadow tone curve
			int Y = CLIP((int)(0.299*r + 0.587*g + 0.114*b));
			tonefactor = (Y>0 ? (float)shtonecurve[Y]/Y : 1);
			r *= tonefactor;
			g *= tonefactor;
			b *= tonefactor;
			
			//brightness/contrast and user tone curve
			r = tonecurve[CLIP(r)];
			g = tonecurve[CLIP(g)];
			b = tonecurve[CLIP(b)];

			if (abs(sat)>0.5 || hCurveEnabled || sCurveEnabled || vCurveEnabled) {
				rgb2hsv(r,g,b,h,s,v);
				if (sat > 0.5) {
					s = (1-(float)sat/100)*s+(float)sat/100*(1-SQR(SQR(1-s)));
				} else {
					if (sat < -0.5)
						s *= 1+(float)sat/100;	
				}
				//HSV equalizer
				if (hCurveEnabled) {
					h = (hCurve->getVal((double)h) - 0.5) * 2 + h;
					if (h > 1.0)
						h -= 1.0;
					else if (h < 0.0)
						h += 1.0;
				}
				if (sCurveEnabled) {
					//shift saturation
					float satparam = (sCurve->getVal((double)h)-0.5) * 2;
					if (satparam > 0.00001) {
						s = (1-satparam)*s+satparam*(1-SQR(1-s));
					} else {
						if (satparam < -0.00001)
							s *= 1+satparam;
					}

					/*s = sCurve->getVal((double)s);
					if (s > 1.0)
						s -= 1.0;
					else if (s < 0.0)
						s += 1.0;*/
				}
				if (vCurveEnabled) {
					//shift value
					float valparam = vCurve->getVal((double)h)-0.5;
					valparam *= (1-SQR(SQR(1-s)));
					if (valparam > 0.00001) {
						v = (1-valparam)*v+valparam*(1-SQR(1-v));
					} else {
						if (valparam < -0.00001)
							v *= (1+valparam);
					}

					/*v = vCurve->getVal((double)v);
					if (v > 1.0)
						v -= 1.0;
					else if (v < 0.0)
						v += 1.0;*/
				}
				hsv2rgb(h,s,v,r,g,b);
			}
			//hsv2rgb(h,s,v,r,g,b);


            int x = (toxyz[0][0] * r + toxyz[1][0] * g + toxyz[2][0] * b) >> 15;
            int y = (toxyz[0][1] * r + toxyz[1][1] * g + toxyz[2][1] * b) >> 15;
            int z = (toxyz[0][2] * r + toxyz[1][2] * g + toxyz[2][2] * b) >> 15;

            x = CLIPTO(x,0,2*65536-1);
            y = CLIPTO(y,0,2*65536-1);
            z = CLIPTO(z,0,2*65536-1);

            int L = cacheL[y];
            lab->L[i][j] = L;
            lab->a[i][j] = CLIPC(((cachea[x] - cachea[y]) * chroma_scale) >> 15);
            lab->b[i][j] = CLIPC(((cacheb[y] - cacheb[z]) * chroma_scale) >> 15);
        }
    }
	
	if (hCurveEnabled) delete hCurve;
	if (sCurveEnabled) delete sCurve;
	if (vCurveEnabled) delete vCurve;
	delete [] cossq;
	//delete [] my_tonecurve;
 }

void ImProcFunctions::luminanceCurve (LabImage* lold, LabImage* lnew, int* curve) {

    int W = lold->W;
    int H = lold->H;
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++)
            lnew->L[i][j] = curve[lold->L[i][j]];
}
		
	
void ImProcFunctions::chrominanceCurve (LabImage* lold, LabImage* lnew, float* acurve, float* bcurve) {
	
	int W = lold->W;
	int H = lold->H;
	/*for (int i=0; i<H; i++)
		for (int j=0; j<W; j++) {
			lnew->a[i][j] = acurve[lold->a[i][j]+32768]-32768;
			lnew->b[i][j] = bcurve[lold->b[i][j]+32768]-32768;
	}*/
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
#pragma omp parallel for if (multiThread)
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++) {
			
			int oa = lold->a[i][j];
			int ob = lold->b[i][j];
			
			float atmp = acurve[oa+32768]-32768;
			float btmp = bcurve[ob+32768]-32768;

            double real_c = 1.0;
            if (params->labCurve.avoidclip) {
				double Lclip = MIN(lnew->L[i][j]/655.35,100.0);
                double cr = tightestroot (Lclip, (double)atmp/chroma_scale, (double)btmp/chroma_scale, 3.079935, -1.5371515, -0.54278342);
                double cg = tightestroot (Lclip, (double)atmp/chroma_scale, (double)btmp/chroma_scale, -0.92123418, 1.87599, 0.04524418);
                double cb = tightestroot (Lclip, (double)atmp/chroma_scale, (double)btmp/chroma_scale, 0.052889682, -0.20404134, 1.15115166);
                if (cr>0 && cr<real_c) real_c = cr;
                if (cg>0 && cg<real_c) real_c = cg;
                if (cb>0 && cb<real_c) real_c = cb;
				//if (i%100==50 && j%100==50) printf ("(i,j)=(%d,%d)  c= %f, rmax= %f \n", i, j, c, real_c);//diagnostic
            }
            
            int nna = (int)((atmp) * real_c );
            int nnb = (int)((btmp) * real_c );
			lnew->a[i][j] = CLIPTO(nna,-32000,32000);
			lnew->b[i][j] = CLIPTO(nnb,-32000,32000);
        }
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
}

#include "cubic.cc"

void ImProcFunctions::colorCurve (LabImage* lold, LabImage* lnew) {

    /*double* cmultiplier = new double [181021];

    double boost_a = (params->colorBoost.amount + 100.0) / 100.0;
    double boost_b = (params->colorBoost.amount + 100.0) / 100.0;

    double c, amul = 1.0, bmul = 1.0;
    if (boost_a > boost_b) {
        c = boost_a;
        if (boost_a > 0)
            bmul = boost_b / boost_a;
    }
    else {
        c = boost_b;
        if (boost_b > 0)
            amul = boost_a / boost_b;
    }

    if (params->colorBoost.enable_saturationlimiter && c>1) {
        // re-generate color multiplier lookup table
        double d = params->colorBoost.saturationlimit * chroma_scale  / 3.0;
        double alpha = 0.5;
        double threshold1 = alpha * d;
        double threshold2 = c*d*(alpha+1.0) - d;
        for (int i=0; i<=181020; i++) { // lookup table stores multipliers with a 0.25 chrominance resolution
            double chrominance = (double)i/4;
            if (chrominance < threshold1)
                cmultiplier[i] = c;
            else if (chrominance < d)
                cmultiplier[i] = (c / (2.0*d*(alpha-1.0)) * (chrominance-d)*(chrominance-d) + c*d/2.0 * (alpha+1.0) ) / chrominance;
            else if (chrominance < threshold2) 
                cmultiplier[i] = (1.0 / (2.0*d*(c*(alpha+1.0)-2.0)) * (chrominance-d)*(chrominance-d) + c*d/2.0 * (alpha+1.0) ) / chrominance;
            else
                cmultiplier[i] = 1.0;
        }
    }
    
	float eps = 0.001;
    double shift_a = params->colorShift.a * chroma_scale + eps, shift_b = params->colorShift.b * chroma_scale + eps;

    short** oa = lold->a;
    short** ob = lold->b;

	#pragma omp parallel for if (multiThread)
    for (int i=0; i<lold->H; i++)
        for (int j=0; j<lold->W; j++) {

            double wanted_c = c;
            if (params->colorBoost.enable_saturationlimiter && c>1) {
                int chroma = (int)(4.0 * sqrt((oa[i][j]+shift_a)*(oa[i][j]+shift_a) + (ob[i][j]+shift_b)*(ob[i][j]+shift_b)));
                wanted_c = cmultiplier [MIN(chroma,181020)];
            }

            double real_c = wanted_c;
            if (wanted_c >= 1.0 && params->colorBoost.avoidclip) {
                double cclip = 100000;
                double cr = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)/chroma_scale*amul, (double)(ob[i][j]+shift_b)/chroma_scale*bmul, 3.079935, -1.5371515, -0.54278342);
                double cg = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)/chroma_scale*amul, (double)(ob[i][j]+shift_b)/chroma_scale*bmul, -0.92123418, 1.87599, 0.04524418);
                double cb = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)/chroma_scale*amul, (double)(ob[i][j]+shift_b)/chroma_scale*bmul, 0.052889682, -0.20404134, 1.15115166);
                if (cr>1.0 && cr<cclip) cclip = cr;
                if (cg>1.0 && cg<cclip) cclip = cg;
                if (cb>1.0 && cb<cclip) cclip = cb;
                if (cclip<100000) {
                    real_c = -cclip + 2.0*cclip / (1.0+exp(-2.0*wanted_c/cclip));
                    if (real_c<1.0)
                        real_c = 1.0;
                }
            }
            
            int nna = (int)((oa[i][j]+shift_a) * real_c * amul);
            int nnb = (int)((ob[i][j]+shift_b) * real_c * bmul);
            lnew->a[i][j] = CLIPTO(nna,-32000,32000);
            lnew->b[i][j] = CLIPTO(nnb,-32000,32000);
        }

    delete [] cmultiplier;*/
}
	
	void ImProcFunctions::impulsedenoise (LabImage* lab) {
		
		if (params->impulseDenoise.enabled && lab->W>=8 && lab->H>=8)
			
			impulse_nr (lab->L, lab->L, lab->W, lab->H, (float)params->impulseDenoise.thresh/20.0 );
	}
	
	void ImProcFunctions::defringe (LabImage* lab) {
		
		if (params->defringe.enabled && lab->W>=8 && lab->H>=8)
			
			PF_correct_RT(lab, lab, params->defringe.radius, params->defringe.threshold, false /*edges only*/ );
	}
	
	void ImProcFunctions::dirpyrdenoise (LabImage* lab) {
		
		if (params->dirpyrDenoise.enabled && lab->W>=8 && lab->H>=8)
			
			dirpyrLab_denoise(lab, lab, params->dirpyrDenoise.luma, params->dirpyrDenoise.chroma, params->dirpyrDenoise.gamma/3.0 );
	}
	
	void ImProcFunctions::dirpyrequalizer (LabImage* lab) {
		
		if (params->dirpyrequalizer.enabled && lab->W>=8 && lab->H>=8) {
			
			//dirpyrLab_equalizer(lab, lab, params->dirpyrequalizer.mult);
			dirpyr_equalizer(lab->L, lab->L, lab->W, lab->H, params->dirpyrequalizer.mult);

		}
	}

void ImProcFunctions::lumadenoise (LabImage* lab, int** b2) {

    /*if (params->lumaDenoise.enabled && lab->W>=8 && lab->H>=8)
#ifdef _OPENMP
#pragma omp parallel
#endif
    	bilateral<unsigned short, unsigned int> (lab->L, lab->L, (unsigned short**)b2, lab->W, lab->H, \
												params->lumaDenoise.radius / scale, params->lumaDenoise.edgetolerance, multiThread);
	 */
}

void ImProcFunctions::colordenoise (LabImage* lab, int** b2) {
	/*
  if (params->colorDenoise.enabled && lab->W>=8 && lab->H>=8) {
#ifdef _OPENMP
#pragma omp parallel
#endif
	  {
	  AlignedBuffer<double>* buffer = new AlignedBuffer<double> (MAX(lab->W,lab->H));
      gaussHorizontal<short> (lab->a, lab->a, buffer, lab->W, lab->H, params->colorDenoise.amount / 10.0 / scale, multiThread);
      gaussHorizontal<short> (lab->b, lab->b, buffer, lab->W, lab->H, params->colorDenoise.amount / 10.0 / scale, multiThread);
      gaussVertical<short>   (lab->a, lab->a, buffer, lab->W, lab->H, params->colorDenoise.amount / 10.0 / scale, multiThread);
      gaussVertical<short>   (lab->b, lab->b, buffer, lab->W, lab->H, params->colorDenoise.amount / 10.0 / scale, multiThread);

      delete buffer;
  }
  }
	*/
}

void ImProcFunctions::getAutoExp  (unsigned int* histogram, int histcompr, double expcomp, double clip, double& br, int& bl) {

    double sum = 0;
    for (int i=0; i<65536>>histcompr; i++)
        sum += histogram[i];

    // compute clipping points based on the original histograms (linear, without exp comp.)
    int clippable = (int)(sum * clip);
    int clipped = 0;
    int aw = (65536>>histcompr) - 1;
    while (aw>1 && histogram[aw]+clipped <= clippable) {
        clipped += histogram[aw];
        aw--;
    }

    clipped = 0;
    int shc = 0;
    while (shc<aw-1 && histogram[shc]+clipped <= clippable) {
        clipped += histogram[shc];
        shc++;
    }

    aw <<= histcompr;
    shc <<= histcompr;
    
    double corr = pow(2.0, expcomp);

    // black point selection is based on the linear result (yielding better visual results)
    bl = (int)(shc * corr);
    // compute the white point of the exp. compensated gamma corrected image
    double awg = (int)(CurveFactory::gamma2 (aw * corr / 65536.0) * 65536.0);

    // compute average intensity of the exp compensated, gamma corrected image
    double gavg = 0;
    for (int i=0; i<65536>>histcompr; i++) 
        gavg += histogram[i] * CurveFactory::gamma2((int)(corr*(i<<histcompr)<65535 ? corr*(i<<histcompr) : 65535)) / sum;

    if (bl < gavg) {
        int maxaw = (gavg - bl) * 4 / 3 + bl; // dont let aw be such large that the histogram average goes above 3/4
        //double mavg = 65536.0 / (awg-bl) * (gavg - bl);
        if (awg < maxaw)
            awg = maxaw;
    }
	
	awg = CurveFactory::igamma2 ((float)(awg/65535.0)) * 65535.0; //need to inverse gamma transform to get correct exposure compensation parameter

	bl = (int)((65535*bl)/awg);
    br = log(65535.0*corr / (awg)) / log(2.0);
    if (br<0) br = 0;
	if (br>10) br=10;
}
	
}

