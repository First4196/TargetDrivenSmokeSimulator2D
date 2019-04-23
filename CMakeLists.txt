project(TargetDrivenSmokeSimulator2D)
 
cmake_minimum_required(VERSION 3.12)

# allows relative library path
cmake_policy(SET CMP0015 NEW)

add_subdirectory(vecmath)

set(SRCS
    fluid_sim.cpp
    linalg.cpp
    main.cpp
)

set(HEADERS
    fluid_sim.hpp
    linalg.hpp
    MAC.h
)

add_executable(TargetDrivenSmokeSimulator2D ${SRCS} ${HEADERS})
target_include_directories(TargetDrivenSmokeSimulator2D PUBLIC GL/include)
target_link_libraries(TargetDrivenSmokeSimulator2D vecmath)
target_link_libraries(TargetDrivenSmokeSimulator2D ${CMAKE_SOURCE_DIR}/GL/lib/Glaux.lib)
target_link_libraries(TargetDrivenSmokeSimulator2D ${CMAKE_SOURCE_DIR}/GL/lib/GLU32.lib)
target_link_libraries(TargetDrivenSmokeSimulator2D ${CMAKE_SOURCE_DIR}/GL/lib/glui32.lib)
target_link_libraries(TargetDrivenSmokeSimulator2D ${CMAKE_SOURCE_DIR}/GL/lib/glut32.lib)
target_link_libraries(TargetDrivenSmokeSimulator2D ${CMAKE_SOURCE_DIR}/GL/lib/OPENGL32.lib)

add_custom_command(TARGET TargetDrivenSmokeSimulator2D POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${PROJECT_SOURCE_DIR}/lib"
    $<TARGET_FILE_DIR:TargetDrivenSmokeSimulator2D>)