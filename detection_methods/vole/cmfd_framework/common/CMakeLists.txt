vole_module_name("common")
vole_module_description("Common classes used in the vole framework")
vole_module_variable("Vole_Common")

vole_add_required_dependencies()
vole_add_optional_dependencies("BOOST_PROGRAM_OPTIONS" "OPENCV" "BOOST_SYSTEM")

vole_compile_library(
    "hashes"
    "vole_config"
    "filesystem_helpers"
    "string_conversion"
    "stopwatch"
    "progress_observer"
)

vole_add_python_module("_common" "pyvole_common.cpp")

vole_add_module()
