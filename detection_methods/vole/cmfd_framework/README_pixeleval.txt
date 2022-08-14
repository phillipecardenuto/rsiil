pixeleval is a small helper program to perform pixelwise evaluation and to
compute common error metrics like false positives, false negatives,... on a
single image.

The code is located in cmfd_framework/pixeleval/

To compile the code, you require
- OpenCV (tested with 2.4.0)

Edit the Makefile such that the header files and libraries of OpenCV are found.
Type "make" to build the code.
Under Windows, setup a Visual Studio solution and configure it properly, such
that OpenCV is included in the build.


The resulting binary "pixeleval" requires a number of parameters in a fixed order:
- ground truth image
- detection result
- the output formatting, denoted as "table" or "human"
- a list of metrics that should be printed/computed, from the set
  {nposground, nnegground, nfp, nfn, ntp, ntn, nposdetected, nnegdetected, tpr, tnr, fpr, fnr, precision, recall, f1, accuracy}

Thus, assume that the ground truth is in an image "red_tower_gt_gj100.png", the
detection results are in a file "red_tower_gj100_matrix.png", and you are
interested in pixelwise precision and recall in a table format, type

./pixeleval red_tower_gt_gj100.png red_tower_gj100_matrix.png table precision recall


The parameters are explained below in greater detail
- {table|human}:
  "human" is a pretty printer, that outputs a list of metrics together with a
  (very brief) description.
  "table" simply puts the results in a table for further processing.

- list of metrics:
  * nposground: The number of positive (tampered) pixels in the ground truth map
  * nnegground: The number of negative (original/untampered) pixels in the
                ground truth map
  * nfp:        false positives, i.e. number of pixels that a) are denoted in
				the ground truth as original and b) are marked by the code as
				tampered.
  * nfn:        false negatives, i.e. number of pixels that a) are denoted in
				the ground truth as tampered and b) are marked by the code as
				original.
  * ntp:        true positives, i.e. number of pixels that a) are denoted in
				the ground truth as tampered and b) are also marked by the code
				as tampered.
  * ntn:        true negatives, i.e. number of pixels that a) are denoted in
				the ground truth as original and b) are also marked by the code
				as original.
  * nposdetected: number of pixels that were marked as positive (tampered) in
                the ground truth map
  * nnegdetected: number of pixels that were marked as negative (original) in
                the ground truth map
  * tpr:        true positive rate, i.e., a value between 0 and 1 relative to
                nposground
  * tnr:        true negative rate, i.e., a value between 0 and 1 relative to
                nnegground
  * fpr:        false positive rate, i.e., a value between 0 and 1 relative to
                nnegground
  * fnr:        false negative rate, i.e., a value between 0 and 1 relative to
                nposground
  * precision:  A popular metric, ntp / (ntp + nfp)
  * recall:     Analogously to precision, ntp / (ntp + nfn)
  * f1:         A joint measure of precision and recall. We used the "plain
                f"-metric, i.e. f1: 2*precision*recall / (precision+recall)
  * accuracy:   (ntp + ntn) / (ntp + nfp + nfn + ntn)

In our performance evaluation of different Copy-Move descriptors, we mainly
relied on Precision, Recall and the F1 measure. However, in order to
incorporate further metrics, the pixeleval code is
directly extensible (search for the word "add", it occurs at four locations in
the code).

