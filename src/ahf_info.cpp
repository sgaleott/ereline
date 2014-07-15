#include "ahf_info.hpp"

std::vector<Od_t>
build_od_list(const std::vector<Pointing_t> & pointings)
{
    std::vector<Od_t> result;
    if(pointings.empty())
	return result;

    Pointing_t first_pointing = pointings.front();
    Od_t current_od { first_pointing.od, first_pointing.id, first_pointing.id, 0 };

    for(auto & pnt : pointings) {
	if(current_od.od == pnt.od) {
	    current_od.last_pointing_id = pnt.id;
	    current_od.num_of_pointings++;
	} else {
	    result.push_back(current_od);

	    current_od.od = pnt.od;
	    current_od.first_pointing_id = pnt.id;
	    current_od.last_pointing_id = pnt.id;
	    current_od.num_of_pointings = 1;
	}
    }

    result.push_back(current_od);
    return result;
}
