#ifndef POST_PROCESS_CONFIG_H
#define POST_PROCESS_CONFIG_H

#include "vole_config.h"

class PPConfig : public vole::Config 
{
public:
	PPConfig(const std::string& prefix = std::string());
	~PPConfig();

	// --- GENERAL
	bool graphical;
	/// cmfd output file (can be the matrix or correlation map, etc)
	std::string inputfile;
	/// the original file at which cmfd has been run
	std::string orig;
	/// additional matrix which contains the matched points (cmfd-output matrix)
	/// normaly only be useful if you pass a correlation map as inputfile
	std::string matrix;
	/// the output directory
	std::string outputdir;

	// --- POST-PROCESSING
	/// remove small areas either by an area threshold and/or by 'useMatches'
	bool rmSmallAreas;
	/// area threshold - duplicated areas must have at least the size:
	/*! if areaThreshold (0,1] areaSize = areaThreshold * image.width * image.height / marks the shift vector
	 *	else areaSize = areaThreshold
	 */
	double areaThreshold;
	/// the morphological operations which should be applied
	std::vector<std::string> morphOp;
	/// how many iterations each morphological operation transformation
	/// should be applied
	std::vector<int> morphIter;
	/*! filles the holes by regarding only the most extreme
	 *	outer contour and filling that one
	 */
	bool fillHolesByContour;

public:
	virtual std::string getString() const;

//    #ifdef VOLE_GUI
//        virtual QWidget *getConfigWidget();
//        virtual void updateValuesFromWidget();
//    #endif// VOLE_GUI

	#ifdef WITH_BOOST
		void printConfig(void);
	#endif // WITH_BOOST	
protected:
    #ifdef WITH_BOOST
		virtual void initBoostOptions();
    #endif // WITH_BOOST

}; 

#endif // COPYMOVEFRAMEWORK_CONFIG_H
