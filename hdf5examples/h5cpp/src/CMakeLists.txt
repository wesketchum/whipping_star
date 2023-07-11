install(FILES ${MY_PUBLIC_HEADERS}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include
    COMPONENT devel)

add_executable              (h5cppwrite h5cppwrite.cpp)

target_link_libraries (h5cppwrite LINK_PUBLIC ${libraries})
install(TARGETS h5cppwrite
    DESTINATION ${CMAKE_INSTALL_PREFIX/bin}
    RUNTIME DESTINATION bin
    )
