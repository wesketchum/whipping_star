include_directories(${CMAKE_SOURCE_DIR}/inc)

add_executable (unit1 unit1.cxx)
target_link_libraries (unit1 SBNfitlib tinyxmllib EventWeight)


file(COPY ${CMAKE_SOURCE_DIR}/unit/unit1a.root DESTINATION ${CMAKE_SOURCE_DIR}/build/unit/)
file(COPY ${CMAKE_SOURCE_DIR}/unit/unit1b.root DESTINATION ${CMAKE_SOURCE_DIR}/build/unit/)
file(COPY ${CMAKE_SOURCE_DIR}/unit/unit1c.root DESTINATION ${CMAKE_SOURCE_DIR}/build/unit/)

file(COPY ${CMAKE_SOURCE_DIR}/unit/unit1a_matrix.root DESTINATION ${CMAKE_SOURCE_DIR}/build/unit/)
file(COPY ${CMAKE_SOURCE_DIR}/unit/unit1b_matrix.root DESTINATION ${CMAKE_SOURCE_DIR}/build/unit/)
file(COPY ${CMAKE_SOURCE_DIR}/unit/unit1c_matrix.root DESTINATION ${CMAKE_SOURCE_DIR}/build/unit/)

file(COPY ${CMAKE_SOURCE_DIR}/xml/unit/unit1a.xml DESTINATION ${CMAKE_SOURCE_DIR}/build/unit/)
file(COPY ${CMAKE_SOURCE_DIR}/xml/unit/unit1b.xml DESTINATION ${CMAKE_SOURCE_DIR}/build/unit/)
file(COPY ${CMAKE_SOURCE_DIR}/xml/unit/unit1c.xml DESTINATION ${CMAKE_SOURCE_DIR}/build/unit/)

