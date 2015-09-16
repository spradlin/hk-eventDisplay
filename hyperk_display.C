{
	gROOT->ProcessLine(".L $WCSIMDIR/libWCSimRoot.so");
	gSystem->AddIncludePath(" -I$FITQUN_ROOT");
	gSystem->AddIncludePath(" -I$WCSIMDIR/include/");
	gROOT->ProcessLine(".L hyperk_esd_html_summary.C+");
	gROOT->ProcessLine(".L THKGamma.C+");
	Long_t* id;
	Long_t* size;
	Long_t* flags;
	Long_t* modtime;
	if(!gSystem->GetPathInfo("../hk-fitqun",id,size,flags,modtime))
		{
			#define FITQUNEXISTS
			gROOT->ProcessLine(".L fitQunDisplay.C+");
		}
		else
		{
			#undef FITQUNEXISTS
		}
			
	gROOT->ProcessLine(".L Picker.C+");
	gROOT->ProcessLine(".x hyperk_esd.C");
}

