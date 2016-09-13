# empty variable to keep track of the `make example` dependencies
set(examples_depends)

add_subdirectory(programs)

# This is the only one requiring ROOT and BAT

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT)
if(ROOT_FOUND)
    find_package(BAT)
    if(BAT_FOUND)
        add_subdirectory(bat_gen)
    endif(BAT_FOUND)
endif(ROOT_FOUND)

# Create custom 'examples' target
add_custom_target(examples DEPENDS ${examples_depends})
