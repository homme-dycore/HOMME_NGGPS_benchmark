# Use make install to copy files to the build directory

file(GLOB ASCII "*.ascii") # get ascii files
install(PROGRAMS ${ASCII} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
