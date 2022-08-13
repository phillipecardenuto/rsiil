#ifndef POST_PROCESS_COMMAND_H
#define POST_PROCESS_COMMAND_H

#include "command.h"
#include "post_process_config.h"
#include <opencv2/core/core.hpp>

class PostProcessCommand : public vole::Command 
{
public:
	PostProcessCommand(void);
	~PostProcessCommand(void) { /*delete i;*/ }
	/// execution with saving
	int execute(void);
	/// execution without saving
	cv::Mat* execute_headless(void);
	
	void printShortHelp(void) const;
	void printHelp(void) const;

	void printConfig(void);
	
	// this is our config struct. please keep it public, baby!
	PPConfig config;
	// temporary image 	 
	cv::Mat *i;
};

#endif // POST_PROCESS_COMMAND_H
