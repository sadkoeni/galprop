macro(CHECK_CXX_NATIVE_REGEX_WORKS NATIVE_REGEX_WORKS_VAR)
    set(source_file_name "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeTmp/regex.cxx")

    file(WRITE ${source_file_name}.in
    "
        #include <regex>
        int main() {
            try {
                std::regex re(\"((?:a[bc]d)|(?:def)|(?:e[fg]h))\");
                std::string s(\"abcd.lmno.defg.file\");
                if (std::regex_search(s, re)) {
                    return 0;
                }
            } 
            catch (...) {
                return 2;
            }
            return 1;
        }
    ")

    configure_file(${source_file_name}.in ${source_file_name})

    try_run(
      run_result
      compile_result
      ${CMAKE_BINARY_DIR}
      ${source_file_name}
      COMPILE_OUTPUT_VARIABLE compile_output
      RUN_OUTPUT_VARIABLE run_output
      )

     message("compile_result='${compile_result}'")
     message("run_result='${run_result}'")
     message("compile_output='${compile_output}'")
     message("run_output='${run_output}'")

    if(NOT ${compile_result}) 
        set(${NATIVE_REGEX_WORKS_VAR} 0)
    elseif("${run_result}" EQUAL 0)
        set(${NATIVE_REGEX_WORKS_VAR} 1)
    else()
        set(${NATIVE_REGEX_WORKS_VAR} 0)
    endif()
endmacro()
