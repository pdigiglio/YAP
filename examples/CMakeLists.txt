# empty variable to keep track of the `make example` dependencies
set(examples_depends)

add_subdirectory(programs)

# This is the only one requiring ROOT and BAT
find_package(BAT)
if(BAT_FOUND)
    find_package(ROOT)
    if(ROOT_FOUND)
        add_subdirectory(bat_gen)
    endif(ROOT_FOUND)
endif(BAT_FOUND)

# Create custom 'examples' target
add_custom_target(examples DEPENDS ${examples_depends})
