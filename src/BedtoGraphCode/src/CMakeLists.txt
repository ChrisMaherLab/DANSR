cmake_minimum_required(VERSION 2.8)

set(LIB_SOURCES
    ALGraph.h
    Alignments.cpp
    Alignments.h
    BedToGraph.cpp
    BedToGraph.h
    ClusterGraph.cpp
    ClusterGraph.h
    Clusters.cpp
    Clusters.h
    ReadAlignments.cpp
    ReadAlignments.h
    AlignmentClusters.cpp
    AlignmentClusters.h
    Gene3.cpp
    Gene3.h
    Utils.cpp
    Utils.h
    Groups.cpp
    Groups.h
    main.cpp
)

add_library(bedtograph ${LIB_SOURCES})
target_link_libraries(bedtograph)
add_executable(bedtograph-bin main.cpp)
target_link_libraries(bedtograph-bin bedtograph)
set_target_properties(bedtograph-bin PROPERTIES OUTPUT_NAME bedToGraph)
