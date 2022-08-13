The ground truth computation consists of C++ code for the actual splicing, and
perl scripts to guide the splicing process.

We recommend to use the perl scripts for convenience. Two variables need to be
adjusted before the scripts can be used: 
- edit the file cmfd_framework/ground_truth_db/scripts/db_setup.pl
- in line 10, enter the full (absolute) path to the vole binary (typically
  something like '<workdir>/cmfd_framework/build/bin/vole')
- in line 12, enter the full (absolute) path to the root directory of the
  dataset (e.g., something like '<workdir>/benchmark_data/')

Ground truth is created with the help of 'configurations'. These configurations
are stored in the file cmfd_framework/ground_truth_db/scripts/splice_configs.pl
A configuration is a key in the hash %splice_config. For instance, 'rot'
creates rotated copies, from 2 degrees to 10 degrees in steps of 2 degrees,
i.e. 6*48=288 output images (for 48 database images).

To create ground truth for the provided dataset, it is most convenient to call
  cmfd_framework/ground_truth_db/scripts/all_calls_splice_image.pl <output_dir> <configuration>
, for instance
  cmfd_framework/ground_truth_db/scripts/all_calls_splice_image.pl /tmp/ rot

If the configuration is omitted (or the script is called without any
parameters), it shows all available configuration entries. Note that the output
directory must already exist.

all_calls_splice_image.pl is a convenience script that repeatedly calls
splice_image.pl. This script again calls then the low-level C++ code.


In case that you require more customized calls, invoke splice_image.pl directly. It expects at least 3 parameters:
-c            is an optional argument, that prints the generated commands to
			  stdout instead of directly executing them. This is either useful
			  for the creation of batch scripts, or for debugging.
<output_dir>  the directory where to store the images. Note that the directory
              must already exist.
<config>      the configuration, according to the configuration list in the
              file splice_configs.pl.
<image_numbers> the image numbers, according to db_setup.pl, that should be
                used. Just add them as separate parameters.

One call, following the example from above, could be
  cmfd_framework/ground_truth_db/scripts/splice_image.pl /tmp/ rot 2 3 4
to create rotated copy-move forgeries for only the images number 2, 3 and 4
(acc. to db_setup.pl, these are 'ship_number', 'tapestry' and 'beach_wood').

To create own ground truth, it is required to create valid dataset entries in
the file db_setup.pl. One dataset entry consists of the keys
- path:                the full absolute path to the image
- orig:                the file name of the original (unmodified) image
- raw_copy:            the file name of a 1-to-1 reference copy, or the output
                       file if copies should be created
- copy:                the file name of the output image
- full_snippets:       a list of copied regions ('snippets') in separate files.
                       The images have the size of the original image.
- full_snippets_alpha: a list of additional black & white images that contain
                       alpha masks for the full_snippets entries. 
- snippets:            file names for the cropped snippet files (i.e., fully
                       transparent lines and columns cropped)
- snippets_alpha:      file names for the alpha channels for the 'snippet'
                       files.
- ground_truth:        file name for the ground truth to this image
- gt_snippets:         file name for the ground truth snippets, used to reverse
                       engineer the snippet position (if insert_positions need to be determined).
- opaque:              threshold for the alpha channel. Intensities above this
					   threshold are considered as fully opaque. In our
					   experiments, we typically set it to 200 to overcome
					   weaknesses in our input data.
- prefix:              the 'important part' of the file name, used to create the output file name.
- insert_positions:    A comma-separated list in the format 'x1,y1,x2,y2,...' of pixel positions
                       where the snippets 1,2,... should be inserted.
- source_positions:    A comma-separated list in the format 'x1,y1,x2,y2,...' of pixel positions
					   where the snippets are taken from. The script
					   'find_snippet_positions.pl' can be used to reverse
					   engineer these positions, in case that it is not clear
					   from the creation process.
- l2dist:              Whether the L2-distance (Euclidean distance) should be
					   used to reverse-engineer the source_positions. If set to
					   0, the L1-distance is used.


Depending on how the ground truth is created, it might be necessary to
reverse-engineer the position where a snippet was extracted. In this case,
find_snippet_position.pl can be used. As an argument, it expects the image
numbers (according to db_setup.pl) of the files to be processed. Output is
written to stdout, and must manually be inserted in the fields
'insert_positions' and 'source_positions' of the entries of db_setup.pl.
