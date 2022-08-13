#include "post_process_config.h"

#ifdef VOLE_GUI
// add qt includes
#endif // VOLE_GI

using namespace boost::program_options;

PPConfig::PPConfig(const std::string & prefix)
	: Config(prefix)
{
#ifdef WITH_BOOST
	initBoostOptions();
#endif // WITH_BOOST
}

PPConfig::~PPConfig() {
}

std::string PPConfig :: getString() const 
{
	std::stringstream s;
	if (prefix_enabled){
		s << "[" << prefix << "]" << std::endl;
	}
	else {
		s << "input=" << inputfile << " # Image to process" << std::endl
		  << "output=" << outputdir << " # Working directory" << std::endl
		  << "orig=" << orig
		  << " # original image in which a possible contour can be drawn" << std::endl
		  << "matrix=" << matrix
		  << " # matrix which contains matched pairs, used to remove wrong areas" << std::endl;
	}
	s << "rmSmallAreas= " << rmSmallAreas
	  << " # remove small areas either by an area threshold and/or by 'useMatches'" << std::endl
	  << "areaThreshold=" << areaThreshold
	  << " # area threshold - duplicated areas must have at least this size" << std::endl
	  << "\tif areaThreshold (0,1] -> size = areaThreshold * image.width * image.height / marks the shift vector\n"
	  << "\telse size = areaThreshold\n"
	  << "fillHolesByContour=" << fillHolesByContour
	  << " # filles the holes by regarding only the most extreme\n"
	  << "\touter contour and filling that one" << std::endl
	  << "morphOp=";
	for (size_t i = 0; i < morphOp.size() ; i++)
		s << morphOp[i] << " ";
	s << " # operation for morphological filter, e.g. ERODE, DILATE, OPEN, CLOSE, etc (see opencv doku)\n";
	s << "morphIter=";
	for (size_t i = 0; i < morphIter.size() ; i++)
		s << morphIter[i] << " ";
	s << " #  iterations for every operation (you can give every operation a different number of iteration)";
	s << std::endl;
	return s.str();
}

#ifdef VOLE_GUI
QWidget * PPConfig :: getConfigWidget() 
{
	// create a qt widget here, fill the widget of the parent class
	this->initConfigWidget();
	QVBoxLayout *cmfd_config = new QVBoxLayout();
	// ...
	layout->addLayout(cmfd_config);
	configWidget->setLayout(layout);
	return configWidget;
}

void PPConfig :: updateValuesFromWidget() 
{
	// pull values out of the widget defined in getConfigWidget()
}
#endif //VOLE_GUI

#ifdef WITH_BOOST
void PPConfig :: initBoostOptions()
{
	if ( ! prefix_enabled ){
		options.add_options()
				(key("graphical"), bool_switch(&graphical)->default_value(false),
				 "Show any graphical output during runtime")
				(key("input,I"), value(&inputfile)->default_value("input.png"),
				 "Image to process")
				(key("orig"), value(&orig)->default_value(""),
				 "original image in which a possible contour can be drawn")
				(key("matrix"), value(&matrix)->default_value(""),
				 "matrix which contains matched pairs, used to remove wrong areas")
				(key("output,O"), value(&outputdir)->default_value("/tmp/"),
				 "Working directory")
				;
	}
	options.add_options()
			(key("rmSmallAreas"), bool_switch(&rmSmallAreas)->default_value(false),
			 "remove small areas either by an area threshold and/or by 'useMatches'")
			(key("areaThreshold"), value(&areaThreshold)->default_value(0.0),
			 "area threshold - duplicated areas must have at least the size:\n"
			 "if areaThreshold (0,1] -> size = areaThreshold * image.width * image.height / marks the shift vector"
			 "else size = areaThreshold")
			(key("fillHolesByContour"), bool_switch(&fillHolesByContour)->default_value(false),
			 "filles the holes by regarding only the most extreme\n"
			 "outer contour and filling that one")
			(key("morphOp"), value(&morphOp),
			 "the morphological operations which should be applied")
			(key("morphIter"), value(&morphIter),
			 "how many iterations each morphological operation transformation\n"
			 "should be applied")
			;
}

void PPConfig :: printConfig(void)
{
    // TODO
}
#endif // WITH_BOOST

