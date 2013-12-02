#ifndef R_HPP
#define R_HPP

#include <pair>
#include "pointing_id.hpp"

typedef float r_type; // Dimensionless number
typedef std::pair<pointing_id, r_type> r_value;

struct r_datatype {
  std::vector<r_value> r_values;
};

#endif
