This directory contains the source code for the copy-move framework, the ground
truth generation and a small program for evaluation, i.e., computing pixelwise
error rates.

The code was the backbone of the paper
"Evaluation of Popular Copy-Move Forgery Detection Approaches" by V.
Christlein, C. Riess, J. Jordan, C. Riess and E. Angelopoulou, IEEE TIFS, vol.
7, no. 6, pp. 1841-1854, Dec. 2012.

-------------------------------------------------------------------------------

The description of the code is split to different readme-files. You might want
to go through it in this order:

README_build.txt         - description of how to build the software
README_ground_truth_examples.txt - description on how to create ground truth images
README_commands.txt      - brief description of the commands in the software package
README_cmfd_examples.txt - example calls for creating ground truth and detecting copy-move forgeries
README_cmfd_options.txt  - more elaborate description of the options for
                           copy-move forgery detection - fine-tuning the results
                           is fairly complex for this command.
README_pixeleval.txt     - description of a small helper program to perform
                           pixelwise performance evaluation


We received helpful comments from various sides.
In particular, we would like to thank
- Li Jian and
- Seung-Jin Ryu
for their bug reports and fixes. Note that the whole code is "research
quality". Nevertheless, comments are welcome, and if time permits, we are also
happy to incorporate fixes, as far as the intended functionality of the code is
concerned (please address them to "Christian dot Riess at fau dot de" or
"Vincent dot Christlein at fau dot de").
