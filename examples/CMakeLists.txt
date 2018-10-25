include(${ROOT_USE_FILE})

include_directories ( ${CMAKE_SOURCE_DIR}/inc)
ROOT_GENERATE_DICTIONARY(G__EventWeight eventweight.h LINKDEF eventweight_LinkDef.h)
#---Create a shared library with geneated dictionary
add_library(EventWeight SHARED eventweight.cxx G__EventWeight.cxx)
target_link_libraries(EventWeight ${ROOT_LIBRARIES})

add_executable (example1 example_1.cxx)
target_link_libraries (example1 SBNfitlib tinyxmllib EventWeight)

add_executable (example2 example_2.cxx)
target_link_libraries (example2 SBNfitlib tinyxmllib EventWeight)

add_executable (example3 example_3.cxx)
target_link_libraries (example3 SBNfitlib tinyxmllib EventWeight)

#add_executable (example4 example_4.cxx)
#target_link_libraries (example4 SBNfitlib tinyxmllib EventWeight)

file(COPY ${CMAKE_SOURCE_DIR}/xml/example.xml DESTINATION ${CMAKE_SOURCE_DIR}/build/examples/)
