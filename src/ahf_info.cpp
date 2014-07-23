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

////////////////////////////////////////////////////////////////////////////////

std::vector<Pointing_t>::const_iterator
find_pid(std::vector<Pointing_t>::const_iterator first,
     std::vector<Pointing_t>::const_iterator last,
     int pid)
{
    decltype(first) it;
    std::iterator_traits<decltype(first)>::difference_type count, step;
    count = std::distance(first, last);

    while(count > 0) {
        it = first;
        step = count / 2;
        std::advance(it, step);
        if(it->id == pid)
            return it;
        else if(it->id < pid) {
            first = ++it;
            count -= step + 1;
        } else count = step;
    }
    return first;
}

////////////////////////////////////////////////////////////////////////////////

void
get_pid_iterators_for_range(const std::vector<Pointing_t> & list_of_pointings,
                            const Range_t<int> & pid_range,
                            std::vector<Pointing_t>::const_iterator & first,
                            std::vector<Pointing_t>::const_iterator & last)
{
    first = find_pid(list_of_pointings.begin(),
                     list_of_pointings.end(),
                     pid_range.start);
    last = find_pid(list_of_pointings.begin(),
                    list_of_pointings.end(),
                    pid_range.end);
}
