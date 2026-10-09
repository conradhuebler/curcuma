# Claude Generated (Oct 2026): helpers to use the structure library (test_cases/structures) from CMake.
# Rules: test_cases/structures/README.md

# curcuma_structure_path(<id> <out_var>): absolute path of the library file of a structure id (xyz or vtf).
function(curcuma_structure_path id out_var)
    file(GLOB _found "${CMAKE_SOURCE_DIR}/test_cases/structures/*/${id}.xyz"
                     "${CMAKE_SOURCE_DIR}/test_cases/structures/*/${id}.vtf")
    list(LENGTH _found _n)
    if(NOT _n EQUAL 1)
        message(FATAL_ERROR "structure '${id}' is not in test_cases/structures exactly once")
    endif()
    set(${out_var} "${_found}" PARENT_SCOPE)
endfunction()

# curcuma_stage_structures(<dest dir> <name>=<id> ...): copy library structures byte for byte into <dest dir> at
# configure time under <name><extension>. Tests read the staged copies, so caches written next to an input
# (for example <name>.topo.json) never end up in the library.
function(curcuma_stage_structures dest)
    file(MAKE_DIRECTORY "${dest}")
    foreach(_item IN LISTS ARGN)
        if(NOT _item MATCHES "^([^=]+)=(.+)$")
            message(FATAL_ERROR "curcuma_stage_structures: expected <name>=<id>, got '${_item}'")
        endif()
        set(_name "${CMAKE_MATCH_1}")
        set(_id "${CMAKE_MATCH_2}")
        curcuma_structure_path("${_id}" _path)
        get_filename_component(_ext "${_path}" EXT)
        string(REGEX REPLACE "^.*(\\.[a-z]+)$" "\\1" _ext "${_ext}")
        configure_file("${_path}" "${dest}/${_name}${_ext}" COPYONLY)
    endforeach()
endfunction()
