/** 
 * \author Johannes Jordan jordan@i5.informatik.uni-erlangen.de
 * \author Vincent Christlein v.christlein@gmail.com
 *
 * ! Extensions / Modifications are welcome !
 *
 * --- Description ---
 *
 * This code may be used for evaluating a binary ground-truth image with 
 * another binary result map to get measures like 'false positive rate', etc.
 *
 * example usage:
 * 1) ./pixeleval groundtruth.png result.png table f1 precision recall
 *	  -> prints out a line with tab-separated 'f1-measure', 'precision' and 'recall'
 *
 * 2) ./pixeleval groundtruth.png result.png human tpr fpr 
 *	  -> prints out in kinda human readable fashion the measures 'true positive
 *			rate' and 'false positive rate'
 */

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>

using namespace cv;

typedef Mat_<float> MAP;

// possible measurements and its human readable presentation
// somewhat ugly, but easy to extend
static const int possible_size = 16;  // edit this if you add more measurements!
static const std::string possible[] = {"nposground", "nnegground", // of ground truth map! 
										"nfp", "nfn", 
										"ntp", "ntn",
										"nposdetected", "nnegdetected", // of result map! 
										"tpr", "tnr", 
										"fpr", "fnr", 
										"precision", "recall", // Note: recall is actually the same as tp-rate
										"f1", "accuracy"
										}; //-> add here more measurements
static const std::string possible_human[] = {	"# total positive (groundtruth)", "# total negative (groundtruth)", 
												"# false positive", "# false negative", 
												"# true positive", "# true negative",
												"# total positive (detected)", "# total negative (detected)",
												"true positive rate", "true negative rate",
												"false positive rate", "false negative rate",
												"precision", "recall",
												"f1-measure", "accuracy"
											}; //-> add here more measurements

// print usage
void usage(std::string program_name){
	std::cerr << "Usage: " << program_name
		 << " <ground_truth_map> <detection_map> <mode> [<mode> = table/human] <measure1...measure2 > [<measureX> =  ";
	for (int i = 0; i < possible_size; i++){
		 std::cerr << possible[i] << ", ";
	}
	std::cerr << "]" << std::endl;
}

int main(int argc, char** argv) {
	// fill possible_which
	if (argc < 5) {
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}
	// read in parameters
	MAP ground = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
	if (ground.empty()) {
		std::cerr << "error reading ground truth image\n";
		exit(EXIT_FAILURE);
	}

	MAP detected = imread(argv[2], CV_LOAD_IMAGE_GRAYSCALE);
	if (detected.empty()) {
		std::cerr << "error reading detection result image\n";
		exit(EXIT_FAILURE);
	}

	std::string mode = argv[3];
	if (mode != "table" && mode != "human"){
		std::cerr << "bad <mode>\n";
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}
	
	// convert and check if argv[4] to argv[argc-1] is a possible string
	// wrong measure-strings will be ignored and raise a warning
	std::vector<std::string> measure;
	for (int i = 4; i < argc; i++){
		bool works = false;
		for (int k = 0; k < possible_size; k++){
			if (argv[i] == possible[k]){
				measure.push_back(argv[i]);
				works = true;
			}
		}
		if (!works)
			std::cerr << "WARNING: " << argv[i] << " will be ignored - not a valid <measure> \n";
	}
	if (measure.empty()){
		std::cerr << "no valid <measures> given\n";
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}
	
	// treat the brightest pixels (typically white) as positive samples
	double mi, ma;
	minMaxLoc(ground, &mi, &ma);
	threshold(ground, ground, ma-1., 1., THRESH_BINARY);

	// treat *any* non-null response of the detector as positive decision
	threshold(detected, detected, 0., 1., THRESH_BINARY);

//	MAP false_positive = max(detected - ground, 0.); does not work OpenCV docs are wrong!
	MAP false_positive = detected - ground;
	max(false_positive, 0., false_positive);

	MAP true_positive = detected - false_positive;
	max(true_positive, 0., true_positive);

	// actual computation
	double total = ground.rows * ground.cols;
	double npositive = sum(ground)[0];
	double nnegative = total - npositive;
	
	double ntruepos = sum(true_positive)[0];
	double nfalsepos = sum(false_positive)[0];
	double ntrueneg = nnegative - nfalsepos;
	double nfalseneg = npositive - ntruepos;
	double nposdetected = ntruepos + nfalsepos;
	double nnegdetected = ntrueneg + nfalseneg;

	double trueposrate = 0;
	if (npositive > 0)
		trueposrate = ntruepos / npositive;
	double truenegrate = 0;
	if (nnegative > 0)
		truenegrate = ntrueneg / nnegative;
	double falseposrate = 0;
	if (nnegative > 0)
		falseposrate = nfalsepos / nnegative;
	double falsenegrate = 0;
	if (npositive > 0)
		falsenegrate = nfalseneg / npositive;
	double precision = 0;
	if ((ntruepos + nfalsepos) > 0)
		precision = ntruepos / (ntruepos + nfalsepos);
	// Note: recall is actually the same as tp-rate
	double recall = 0;
	if ((ntruepos + nfalseneg) > 0)
		recall = ntruepos / (ntruepos + nfalseneg);
	double f1 = 0;
	if ((precision + recall ) > 0)
		f1 = 2*precision*recall / (precision + recall);
	double accuracy = 0;
	if ((ntruepos + nfalseneg + ntrueneg + nfalsepos) > 0)
		accuracy = (ntruepos + ntrueneg) / (ntruepos + nfalseneg + ntrueneg + nfalsepos);
	//-> compute more measures here

	// assign vals to the short-string
	// (unfortunately only boost::assign would make this shorter)
	std::map<std::string, double> val; 
	val[possible[0]] = npositive;
	val[possible[1]] = nnegative;
	val[possible[2]] = nfalsepos;
	val[possible[3]] = nfalseneg;
	val[possible[4]] = ntruepos;
	val[possible[5]] = ntrueneg;
	val[possible[6]] = nposdetected;
	val[possible[7]] = nnegdetected;
	val[possible[8]] = trueposrate;
	val[possible[9]] = truenegrate;
	val[possible[10]] = falseposrate;
	val[possible[11]] = falsenegrate;
	val[possible[12]] = precision;
	val[possible[13]] = recall;
	val[possible[14]] = f1;
	val[possible[15]] = accuracy;
	//-> add more measures to this map here

	// build up a 2nd map of possible-strings to human-readable strings
	std::map<std::string, std::string> human_map;
	for (int i = 0; i < possible_size; i++){
		human_map[possible[i]] = possible_human[i];
	}
	// print results
	for (size_t i = 0; i < measure.size(); i++){
		if (mode == "human"){
			std::cout << human_map[measure[i]] << ": ";
			// don't give it as percent if its actual a number
			if (measure[i].compare(0, 1, "n") == 0)
				std::cout << val[measure[i]] << std::endl;
			else 
				std::cout << val[measure[i]] * 100 << " %" << std::endl;
		}
		else{
			printf("%.4f\t",val[measure[i]]);
		}
	}
	std::cout << std::endl;
}

