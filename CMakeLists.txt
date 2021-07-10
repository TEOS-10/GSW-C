cmake_minimum_required(VERSION 3.12)
project(GSW-C)

add_library(libgswteos-10 SHARED gsw_oceanographic_toolbox.c gsw_saar.c)
add_library(libgswteos-10-static STATIC gsw_oceanographic_toolbox.c gsw_saar.c)
add_executable(gsw_check gsw_check_functions.c gsw_oceanographic_toolbox.c gsw_saar.c)

# Clean data files
if (WIN32)
    add_custom_command(TARGET gsw_check 
                        PRE_BUILD 
                        COMMAND Del gsw_check_data.c
                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/)
    add_custom_command(TARGET gsw_check 
                        PRE_BUILD 
                        COMMAND Del gsw_saar_data.c
                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/)
endif (WIN32)

if (UNIX)
    add_custom_command(TARGET gsw_check 
                        PRE_BUILD 
                        COMMAND rm gsw_check_data.c
                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/)
    add_custom_command(TARGET gsw_check 
                        PRE_BUILD 
                        COMMAND rm gsw_saar_data.c
                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/)
endif (UNIX)

# Create new data
add_custom_command(TARGET gsw_check 
                    PRE_BUILD 
                    COMMAND python ${PROJECT_SOURCE_DIR}/make_check_data.py
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/)
add_custom_command(TARGET gsw_check 
                    PRE_BUILD 
                    COMMAND python ${PROJECT_SOURCE_DIR}/make_saar_data.py
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/)
