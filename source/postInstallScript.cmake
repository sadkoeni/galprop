#Get the link program
execute_process(COMMAND "${CMAKE_COMMAND} -E create_symlink galprop gpavXCO" WORKING_DIR ${CMAKE_INSTALL_PREFIX}/${INSTALL_BIN_DIR})
