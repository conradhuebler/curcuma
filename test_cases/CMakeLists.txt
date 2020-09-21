file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input_aa.xyz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input_ab.xyz
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(rmsd_test
        rmsd/main.cpp)

target_link_libraries(rmsd_test curcuma_core)

add_executable(reorder_test
        reorder/main.cpp)

target_link_libraries(reorder_test curcuma_core)
