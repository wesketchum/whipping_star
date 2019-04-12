include_directories(${CMAKE_SOURCE_DIR}/inc)

add_executable (example1 example_1.cxx)
target_link_libraries (example1 SBNfitlib tinyxmllib EventWeight)

add_executable (example2 example_2.cxx)
target_link_libraries (example2 SBNfitlib tinyxmllib EventWeight)

add_executable (example3 example_3.cxx)
target_link_libraries (example3 SBNfitlib tinyxmllib EventWeight)

add_executable (exampletest example_test.cxx)
target_link_libraries (exampletest SBNfitlib tinyxmllib EventWeight)


#add_executable (exampleI example_i.cxx)
#target_link_libraries (exampleI SBNfitlib tinyxmllib EventWeight)

#add_executable (exampleI2 example_i2.cxx)
#target_link_libraries (exampleI2 SBNfitlib tinyxmllib EventWeight)

#add_executable (example4 example_4.cxx)
#target_link_libraries (example4 SBNfitlib tinyxmllib EventWeight)

file(COPY ${CMAKE_SOURCE_DIR}/xml/example.xml DESTINATION ${CMAKE_SOURCE_DIR}/build/examples/)



add_subdirectory(MiniBooNE_suite)

