#ifndef TSHARCANALYSIS_H
#define TSHARCANALYSIS_H

#include <TObject.h>
#include <TMath.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TList.h>
#include <TCanvas.h>
#include <TH2F.h>

#include "TSRIM.h"
#include "TSharc.h"
#include "TReaction.h"

#ifndef PI
#define PI                       (TMath::Pi())
#endif

#ifndef R2D
#define R2D                       (TMath::RadToDeg())
#endif

#ifndef D2R
#define D2R                       (TMath::DegToRad())
#endif

#ifndef SHCDIR
#define SHCDIR "/Users/steffencruz/Desktop/Steffen/Work/PhD/TRIUMF/CodesAndTools/SharcAnalysis"
#endif

class TSharcAnalysis 	{
	
	public:

		TSharcAnalysis(double Xoff = 0.00, double Yoff = 0.00, double Zoff = 0.00);
		~TSharcAnalysis();
	
		static void Print(Option_t * = "");			
		static void Clear(Option_t * = "");		
		static void ResetDetectors();

		static TVector3 GetPosition(int det, int fs, int bs) 
					{return TSharc::GetPosition(det,fs,bs)+position_offset;}

		static void SetTarget(double x=0.0, double y=0.0, double z=0.0, double thickness=0.0, const char *material="",double dx=0.0, double dy=0.0, double dz=0.0);
		static void SetFracTarg(double frac) {frac_targ = frac;}
		static void SetTargetPosition(TVector3 vec){ position_offset = vec; }		
		static void SetTargetPositionError(TVector3 vec){ position_error = vec; }		
		
    static void SetDelThickness(int det, double t) { DELthicknesses[det-1] = t; }
    static void SetPadThickness(int det, double t) { PADthicknesses[det-1] = t; }
    static void SetDelDeadlayer(int det, double t) { DELdeadlayers[det-1] = t; }
    static void SetPadDeadlayer(int det, double t) { PADdeadlayers[det-1] = t; }
        
    // pixel corners and solid angles in lab frame
    static TVector3 GetEdgePoints(int edgenum, int det, int fs, int bs); //!  
    static double GetSumSolidAngle(int det_min, int fs_min, int bs_min, int det_max = 0, int fs_max = -1, int bs_max = -1); //!
    static double GetSolidAngle(int det, int fs, int bs); //!
    static TH2F *GetSolidAngleMatrix(int det, Bool_t badstrips=false); //!
    
    static void SetOmegaMaxMin(bool printout = false); //!
    static double GetOmegaMax() { return Omega_max; }; //!
    static double GetOmegaMin() { return Omega_min; }; //!   
            
    // Full acceptance curves and corrections
		static TH1D *SetAcceptance(int nmax=10,TReaction *r=NULL,const char *stripsfilename="");
    static TList *GetAcceptanceList(TReaction *r=NULL,const char *stripsfile="",Int_t resolution=10); //! 
    static double RandomizeThetaLab(int det, int fs, int bs);//! 
    static TH1D *RebinAcceptance(TH1D *h, int binsz, Bool_t verbose=false); //! 
    static TH1D *MakeSin(Double_t binsz=1.0); //!
            
    static TH1D *SimulateAngDist(const char *sname, double &cntres, double &reserr, double t1=0.0, double t2=20.0, double counts=1e4, double cnterr=0); //!
            
    // number of strips
    static int GetFrontStrips(int det){ 	return frontstrips[det-1];		}//!
    static int GetBackStrips(int det) {	  return backstrips[det-1];			}//!
           
    // bad strips
    static void SetBadStrips(std::string stripsfile="");    
    static int BadStrip(int det=-1, int fs=-1, int bs=-1);//! 

    // thicknesses of detectors
		static double GetDetectorThickness(double dist, int det, int fs, int bs)
				{ return GetDetectorThickness(dist,det,GetPosition(det,fs,bs).Theta(),GetPosition(det,fs,bs).Phi()); }	//!			
		static double GetDetectorThickness(double dist, int det, double theta, double phi); //!
		
		static double   GetDelDeadLayer(int det, int fs, int bs)		  
		  	{return GetDetectorThickness(DELdeadlayers[det-1],det,fs,bs);} //!
		static double   GetPadDeadLayer(int det, int fs, int bs)
		 	  {return GetDetectorThickness(PADdeadlayers[det-1],det,fs,bs);} //!
		static double   GetDelThickness(int det, int fs, int bs)	  
		  	{return GetDetectorThickness(DELthicknesses[det-1],det,fs,bs);} //!		  
		static double   GetPadThickness(int det, int fs, int bs)	  
		  	{return GetDetectorThickness(PADthicknesses[det-1],det,fs,bs);} //!

		static double   GetDelDeadLayer(int det, double theta = PI/2.0, double phi = 0.0)		     
	  		{return GetDetectorThickness(DELdeadlayers[det-1],det,theta,phi);}	//!
		static double   GetPadDeadLayer(int det, double theta = PI/2.0, double phi = 0.0)	     
	  		{return GetDetectorThickness(PADdeadlayers[det-1],det,theta,phi);}	//!
		static double   GetDelThickness(int det, double theta = PI/2.0, double phi = 0.0)		     
		 	 {return GetDetectorThickness(DELthicknesses[det-1],det,theta,phi);}	//!
		static double   GetPadThickness(int det, double theta = PI/2.0, double phi = 0.0)	     
		 	 {return GetDetectorThickness(PADthicknesses[det-1],det,theta,phi);}	//!
		
		// target thickness
		static double GetTargetThickness(int det, int fs, int bs, double frac_depth = 0.5)
						{ return GetTargetThickness(GetPosition(det,fs,bs).Theta(),GetPosition(det,fs,bs).Phi(),frac_depth); }	//!			
		static double GetTargetThickness(double theta=0.0, double phi = 0.0, double frac_depth = 0.5);	//!		
		
		// target properties
    static TVector3 GetTargetPosition()    { return position_offset; } //!
    static TVector3 GetTargetPosError()    { return position_error ; } //!
    static double GetTargetDepth()  		   { return targetthickness; } //!       
		static const char *GetTargetMaterial() { return targmat.c_str(); } //!    
    
    // target energy loss
    static double GetBeamEnergyInTarget(const char *beam = "sr95", double ebeam=524e3, double frac_depth_end=0.5); 
		static double GetResEnergyAfterTarget(const char *res, double eres, double frac_depth_beg=0.5);
		
		// takes edel and epad and returns ekin
    static double GetReconstructedEnergy(TVector3 pos, int det, double edel, double epad=0.0, char ion='p');//!
    static double GetReconstructedEnergy(int det, int fs, int bs, double edel, double epad=0.0, char ion='p')
    			{ return GetReconstructedEnergy(GetPosition(det,fs,bs),det,edel,epad,ion); } //!
    			
		// takes ekin and returns edel and epad (NB detector thicknesses must be correct!)   
    static std::vector<double> GetMeasuredEnergy(TVector3 pos, int det, double ekin, char ion='p', std::string opt="", double edel=-1.0);//!	
    static std::vector<double> GetMeasuredEnergy(int det, int fs, int bs, double ekin, char ion='p', std::string opt="", double edel=-1.0)
    			{return GetMeasuredEnergy(GetPosition(det,fs,bs),det,ekin,ion,opt,edel); } //!		 
    
    static TList *SimulateMeasurement(TReaction *r, Int_t resolution=0, Bool_t use_badstrips=false);
    	
		static void InitializeSRIMInputs();   //!
    static TSRIM *GetSRIM(char ion, std::string material); //!
    static TReaction *GetReaction(void) { return reaction; } //!

	private: 

    static TH1D *MakeScaleHist(Double_t binsz, Double_t scale, Double_t err, const char *name);
		static double ThetaCoverage(Double_t *x, Double_t *p);//!
	  static void GetRid(const char *name, Bool_t delete_all = true); // removes object
	
		static TVector3 position_offset; //!
		static TVector3 position_error; //!
    static double targetthickness;   // microns
    static double targetradius;      // microns
    static std::string targmat;			 // material name
    static double frac_targ;        // reaction position in target [fraction]

    static std::string badstripsfile;	// file containing bad strips
    static Bool_t resetbadstrips;
    static UInt_t nbadstrips;
    static Int_t accres;            // resolution of acceptance
    static TList *coveragelist;     // list containing acceptance results
    
		// detector segmentation and pitches
    static int number_of_detectors;
    static int frontstrips[16];
    static int backstrips[16];
    static double frontpitches[16];
    static double backpitches[16];            
		// various thicknesses in um
		static double DELdeadlayers[16]; 
		static double PADdeadlayers[16]; 
		static double DELthicknesses[16]; 
		static double PADthicknesses[16]; 

    // max min solid angle (pixels)
    static double Omega_min;
    static double Omega_max;    

		static TSRIM *p_in_targ;  //!
		static TSRIM *p_in_si;   //!
		static TSRIM *d_in_targ; //!
		static TSRIM *d_in_si;  //!
		static TSRIM *c_in_targ;   //!
		static TSRIM *c_in_si;   //!
		static TSRIM *t_in_targ;   //!
		static TSRIM *t_in_si;   //!		
		static TSRIM *a_in_targ;  //! 
		static TSRIM *a_in_si;   //!
		static TSRIM *beam_in_targ;  //!
		
		static TReaction *reaction;
		        
	ClassDef(TSharcAnalysis,0)
};


#endif
