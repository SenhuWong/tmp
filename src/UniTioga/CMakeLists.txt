

set(SRC ${SRC}
${CMAKE_CURRENT_SOURCE_DIR}/ADT.cpp
${CMAKE_CURRENT_SOURCE_DIR}/test_list.cpp
${CMAKE_CURRENT_SOURCE_DIR}/BlockInfo.cpp
${CMAKE_CURRENT_SOURCE_DIR}/bookkeeping.cpp
${CMAKE_CURRENT_SOURCE_DIR}/CellBrick.cpp
${CMAKE_CURRENT_SOURCE_DIR}/checkContainment.cpp

${CMAKE_CURRENT_SOURCE_DIR}/exchangeboxes.cpp
${CMAKE_CURRENT_SOURCE_DIR}/exchangeDonor.cpp
${CMAKE_CURRENT_SOURCE_DIR}/exchangeSearchData.cpp
${CMAKE_CURRENT_SOURCE_DIR}/holeMap.cpp 
${CMAKE_CURRENT_SOURCE_DIR}/linkList.cpp
${CMAKE_CURRENT_SOURCE_DIR}/MeshBlock.cpp
# ${CMAKE_CURRENT_SOURCE_DIR}/mul_test.cpp
${CMAKE_CURRENT_SOURCE_DIR}/parallelComm.cpp
${CMAKE_CURRENT_SOURCE_DIR}/recursiveADT.cpp
${CMAKE_CURRENT_SOURCE_DIR}/search.cpp
${CMAKE_CURRENT_SOURCE_DIR}/tioga.cpp
${CMAKE_CURRENT_SOURCE_DIR}/TiogaFeeder.cpp
${CMAKE_CURRENT_SOURCE_DIR}/utils.cpp
${CMAKE_CURRENT_SOURCE_DIR}/write.cpp

#These are under development but for symbols let's link them
${CMAKE_CURRENT_SOURCE_DIR}/amr_index.cpp
${CMAKE_CURRENT_SOURCE_DIR}/exchangeDonorAMR.cpp
${CMAKE_CURRENT_SOURCE_DIR}/CartBlock.cpp 
${CMAKE_CURRENT_SOURCE_DIR}/CartGrid.cpp 
${CMAKE_CURRENT_SOURCE_DIR}/CartOpps.cpp

PARENT_SCOPE
)