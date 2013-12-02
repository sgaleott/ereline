#ifndef CALIBRATION_CONSTANT_HPP
#define CALIBRATION_CONSTANT_HPP

#include "pointing_id.hpp"

typedef float cal_constant_type; // Units are K/V
typedef float cal_offset_type; // Units are K

struct cal_constant_value {
  pointing_id pid;
  cal_constant_type inv_gain;
  cal_offset_type offset;
};

struct cal_constant_datatype {
  std::vector<cal_constant_value> cal_constant_values;
};

#endif
