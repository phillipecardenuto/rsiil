vole_module_name("cmfdgt")
vole_module_description("Ground Truth Manager for Copy-Move Forgery Detection")
vole_module_variable("Vole_CMFD_Ground_Truth")

vole_add_required_dependencies("OPENCV" "BOOST" "BOOST_PROGRAM_OPTIONS")
vole_add_required_modules()

vole_compile_library(
	"command_ground_truth.h"             "command_ground_truth.cpp"
	"command_splice.h"                   "command_splice.cpp"
	"command_reverse_engineer_snippet.h" "command_reverse_engineer_snippet.cpp"
	"command_splice_random_block.h"      "command_splice_random_block.cpp"
)

vole_add_command("gt"             "command_ground_truth.h"             "vole::GroundTruth")
vole_add_command("gt_get_snippet" "command_reverse_engineer_snippet.h" "vole::cmfdgt::CommandReverseEngineerSnippet")
vole_add_command("splice"         "command_splice.h"                   "vole::cmfd_gt::CommandSplice")
vole_add_command("gt_postproc"    "command_gt_post_processing.h"       "vole::cmfdgt::CommandGtPostProcessing")
vole_add_command("splice_rb"	  "command_splice_random_block.h"      "vole::cmfd_gt::CommandSpliceRandomBlock")

vole_compile_library(
	forgery
	good_match
	ground_truth_common
	ground_truth_config
	subimage
	alpha_channel
	reverse_engineer_snippet_config
	gt_post_processing_config
	command_gt_post_processing
	splice_random_block_config
	
	splice_config
	splice_core
)

vole_add_module()
