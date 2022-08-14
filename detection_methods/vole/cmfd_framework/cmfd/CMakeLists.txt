vole_module_name("cmfd")
vole_module_description("Copy-Move Forgery Detection by Vince")
vole_module_variable("Vole_CMFD")

vole_add_required_dependencies("OPENCV" "BOOST_PROGRAM_OPTIONS" "BOOST_THREAD" "BOOST_SERIALIZATION" "BOOST_SYSTEM")
vole_add_required_modules()

vole_add_command("cmfd" "copymoveframework.h" "CopyMoveFramework")
vole_add_command("pp" "post_process/post_process_command.h" "PostProcessCommand")

vole_compile_library(
cmfd_util
blockhandling
preproc
wavelet
feature
luo
blur
dct
pctfeat
hu
svd
dwtfeat
fmt
lin
bravo
cps
circlefeat
zernike
kpca
featfactory 
pca
matching 
verification 
fastsats
buildcluster
mark 
copymoveframework_config
log
execution 
copymoveframework
post_process/post_process_config
post_process/post_process_core
post_process/post_process_command
)

vole_add_module()
