#ifndef TYPES_H
#define TYPES_H

#include <cstdint>
#include <iostream>
#include <limits>

typedef unsigned int ui;
typedef double decimal;

typedef uint32_t VertexID;
typedef short NodeID;

const VertexID INVALID_VID(std::numeric_limits<VertexID>::max());
const NodeID INVALID_MID(std::numeric_limits<NodeID>::max());


#endif