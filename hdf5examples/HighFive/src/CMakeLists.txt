install(FILES ${MY_PUBLIC_HEADERS}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include
    COMPONENT devel)

add_executable              (highfivewrite highfivewrite.cpp)

target_link_libraries (highfivewrite LINK_PUBLIC ${libraries})
install(TARGETS highfivewrite
    DESTINATION ${CMAKE_INSTALL_PREFIX/bin}
    RUNTIME DESTINATION bin
    )
