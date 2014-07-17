/*
 * Squeezer - compress LFI detector pointings and differenced data
 * Copyright (C) 2013 Maurizio Tomasi (Planck collaboration)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */

#include "squeezer.hpp"
#include "logging.hpp"
#include "datatypes.hpp"

#include <algorithm>
#include <memory>
#include <cstring>
#include <fstream>

#include <gsl/gsl_poly.h>

//////////////////////////////////////////////////////////////////////

uint8_t
read_uint8(std::istream & in)
{
    char result;
    in.read(&result, sizeof(result));
    if(in.bad() || in.fail())
	throw std::runtime_error("unexpected end of file");

    return (uint8_t) result;
}

//////////////////////////////////////////////////////////////////////

uint16_t
read_uint16(std::istream & in)
{
    uint8_t byte1 = read_uint8(in);
    uint8_t byte2 = read_uint8(in);

    return (((uint16_t) byte1) << 8) + byte2;
}

//////////////////////////////////////////////////////////////////////

uint32_t
read_uint32(std::istream & in)
{
    uint16_t word1 = read_uint16(in);
    uint16_t word2 = read_uint16(in);

    return (((uint32_t) word1) << 16) + word2;
}

//////////////////////////////////////////////////////////////////////

uint64_t
read_uint64(std::istream & in)
{
    uint32_t double_word1 = read_uint32(in);
    uint32_t double_word2 = read_uint32(in);

    return (((uint64_t) double_word1) << 32) + double_word2;
}

//////////////////////////////////////////////////////////////////////

double
read_double(std::istream & in)
{
    double value;
    in.read(reinterpret_cast<char *>(&value), sizeof(value));
    if(in.bad() || in.fail()) {
	throw std::runtime_error("unexpected end of file");
    }

    return value;
}

////////////////////////////////////////////////////////////////////////////////

enum Squeezer_file_type_t { 
    SQZ_NO_DATA,
    SQZ_DETECTOR_POINTINGS, 
    SQZ_DIFFERENCED_DATA 
};

////////////////////////////////////////////////////////////////////////////////

enum Chunk_type_t { 
    CHUNK_DELTA_OBT = 10,
    CHUNK_SCET_ERROR = 11,
    CHUNK_THETA = 12,
    CHUNK_PHI = 13,
    CHUNK_PSI = 14,
    CHUNK_DIFFERENCED_DATA = 15,
    CHUNK_QUALITY_FLAGS = 16
};

////////////////////////////////////////////////////////////////////////////////

struct Squeezer_file_header_t {
    uint8_t file_type_mark[4];
    double floating_point_check;

    uint16_t date_year;
    uint8_t date_month;
    uint8_t date_day;

    uint8_t time_hour;
    uint8_t time_minute;
    uint8_t time_second;

    uint8_t horn;
    uint8_t arm;

    uint16_t od;

    // Since these values are needed to decompress the SCET chunk, it
    // would be better to gather them into a dedicated sub-structure.
    double first_obt;
    double last_obt;
    double first_scet_in_ms;
    double last_scet_in_ms;

    uint32_t number_of_chunks;

    Squeezer_file_header_t(std::istream & in);

    Squeezer_file_type_t get_type() const {
	if(std::strcmp((char *) file_type_mark, "PDP") == 0)
	    return SQZ_DETECTOR_POINTINGS;
	else if(std::strcmp((char *) file_type_mark, "PDD") == 0)
	    return SQZ_DIFFERENCED_DATA;
	else
	    return SQZ_NO_DATA;
    }

    bool is_valid() const {
	if((strcmp((char *) file_type_mark, "PDP") != 0 && 
	    strcmp((char *) file_type_mark, "PDD") != 0) ||
	   date_year < 2013 ||
	   (date_month < 1 || date_month > 12) ||
	   (date_day < 1 || date_day > 31) ||
	   time_hour > 23 ||
	   time_minute > 59 ||
	   time_second > 59 ||
	   floating_point_check != 231250.0 ||
	   (horn < 18 || horn > 28) ||
	   (arm != 0 && arm != 1) ||
	   first_obt >= last_obt ||
	   first_scet_in_ms >= last_scet_in_ms ||
	   number_of_chunks == 0)
	    return false;
	else
	    return true;
    }
};

//////////////////////////////////////////////////////////////////////

Squeezer_file_header_t::Squeezer_file_header_t(std::istream & in)
{
    file_type_mark[0] = read_uint8(in);
    file_type_mark[1] = read_uint8(in);
    file_type_mark[2] = read_uint8(in);
    file_type_mark[3] = read_uint8(in);

    floating_point_check = read_double(in);

    date_year = read_uint16(in);
    date_month = read_uint8(in);
    date_day = read_uint8(in);

    time_hour = read_uint8(in);
    time_minute = read_uint8(in);
    time_second = read_uint8(in);

    horn = read_uint8(in);
    arm = read_uint8(in);
    od = read_uint16(in);
    first_obt = read_double(in);
    last_obt = read_double(in);
    first_scet_in_ms = read_double(in);
    last_scet_in_ms = read_double(in);
    number_of_chunks = read_uint32(in);
}

//////////////////////////////////////////////////////////////////////

struct Error_t {
    double min_abs_error;
    double max_abs_error;
    double mean_abs_error;
    double mean_error;

    Error_t() {
	min_abs_error = 0.0;
	max_abs_error = 0.0;
	mean_error = 0.0;
	mean_abs_error = 0.0;
    }

    void read_from_file(std::istream & in) {
	min_abs_error = read_double(in);
	max_abs_error = read_double(in);
	mean_abs_error = read_double(in);
	mean_error = read_double(in);
    }

    bool is_valid() const {
	return (min_abs_error >= 0.0) && 
	    (mean_abs_error >= 0.0) &&
	    (min_abs_error <= max_abs_error);
    }
};

//////////////////////////////////////////////////////////////////////

struct Squeezer_chunk_header_t {
    uint8_t chunk_mark[4];
    uint64_t number_of_bytes;
    uint32_t number_of_samples;

    uint32_t chunk_type;

    Error_t compression_error;

    Squeezer_chunk_header_t(std::istream & in) {
	chunk_mark[0] = read_uint8(in);
	chunk_mark[1] = read_uint8(in);
	chunk_mark[2] = read_uint8(in);
	chunk_mark[3] = read_uint8(in);

	number_of_bytes = read_uint64(in);
	number_of_samples = read_uint32(in);

	chunk_type = read_uint32(in);

	compression_error.read_from_file(in);
    }

    bool is_valid() const {
	if(chunk_mark[0] != 'C' ||
	   chunk_mark[1] != 'N' ||
	   chunk_mark[2] != 'K' ||
	   chunk_mark[3] != 0 ||
	   number_of_bytes == 0 ||
	   number_of_samples == 0 ||
	   chunk_type < CHUNK_DELTA_OBT || chunk_type > CHUNK_QUALITY_FLAGS)
	    return false;

	return true;
    }
};

////////////////////////////////////////////////////////////////////////////////

class Byte_buffer_t {
public:
    std::vector<uint8_t> buffer;
    size_t cur_position;

    Byte_buffer_t() 
	: buffer(),
	  cur_position(0) {}

    Byte_buffer_t(size_t length, const uint8_t * raw_buffer);

    size_t size() const {
	return buffer.size();
    }

    // Number of bytes that can still be read
    size_t items_left() const {
	return buffer.size() - cur_position;
    }

    uint8_t read_uint8() {
	return buffer.at(cur_position++);
    }

    uint16_t read_uint16();
    uint32_t read_uint32();
    uint64_t read_uint64();
    float read_float();
    double read_double();
};

////////////////////////////////////////////////////////////////////////////////

Byte_buffer_t::Byte_buffer_t(size_t length, const uint8_t * raw_buffer)
{
    buffer.resize(length);
    std::copy(raw_buffer, raw_buffer + length,
	      buffer.begin());

    cur_position = 0;
}

////////////////////////////////////////////////////////////////////////////////

uint16_t
Byte_buffer_t::read_uint16()
{
    uint8_t byte1 = read_uint8();
    uint8_t byte2 = read_uint8();
    return (((uint16_t) byte1) << 8) + byte2;
}

////////////////////////////////////////////////////////////////////////////////

uint32_t
Byte_buffer_t::read_uint32()
{
    uint16_t word1 = read_uint16();
    uint16_t word2 = read_uint16();
    return (((uint32_t) word1) << 16) + word2;
}

////////////////////////////////////////////////////////////////////////////////

uint64_t
Byte_buffer_t::read_uint64()
{
    uint32_t dword1 = read_uint32();
    uint32_t dword2 = read_uint32();
    return (((uint64_t) dword1) << 32) + dword2;
}

////////////////////////////////////////////////////////////////////////////////

float
Byte_buffer_t::read_float()
{
    static_assert(sizeof(float) == 4, 
		  "This code assumes that single-precision floating "
		  "point values are 4 bytes wide.");

    uint32_t value = read_uint32();
    return *(reinterpret_cast<float *>(&value));
}

////////////////////////////////////////////////////////////////////////////////

double
Byte_buffer_t::read_double()
{
    static_assert(sizeof(double) == 8, 
		  "This code assumes that double-precision floating "
		  "point values are 8 bytes wide.");

    uint64_t value = read_uint64();
    return *(reinterpret_cast<double *>(&value));
}

////////////////////////////////////////////////////////////////////////////////

void
rle_decompression(Byte_buffer_t & input_stream,
		  size_t output_size,
		  std::vector<uint32_t> & output)
{
    output.resize(output_size);
    size_t output_idx = 0;
    while(output_idx < output_size) {
	uint32_t count = input_stream.read_uint32();
	uint32_t value = input_stream.read_uint32();

	std::fill(output.begin() + output_idx,
		  output.begin() + output_idx + count,
		  value);

	output_idx += count;
    }
}

////////////////////////////////////////////////////////////////////////////////

struct Frame_t {
    uint8_t num_of_elements;
    std::vector<double> parameters;

    Frame_t()
	: num_of_elements(0),
	  parameters() {}

    Frame_t(uint8_t a_num_of_elements,
	    std::vector<double> a_parameters)
	: num_of_elements(a_num_of_elements),
	  parameters(a_parameters) {}

    void read_from_buffer(Byte_buffer_t & input_buffer);

    bool is_encoded_as_a_polynomial() const {
	return num_of_elements > parameters.size();
    }
};

typedef std::vector<Frame_t> Vector_of_frames_t;

//////////////////////////////////////////////////////////////////////

void
Frame_t::read_from_buffer(Byte_buffer_t & input_buffer)
{
    num_of_elements = input_buffer.read_uint8();

    size_t num_of_parameters = input_buffer.read_uint8();
    parameters.resize(num_of_parameters);
    for(size_t idx = 0; idx < num_of_parameters; ++idx) {
	parameters[idx] = input_buffer.read_float();
    }
}

////////////////////////////////////////////////////////////////////////////////

void
poly_fit_decode(size_t num_of_elements_to_decode,
		Byte_buffer_t & input_buffer,
		std::vector<double> & values)
{
    values.resize(num_of_elements_to_decode);

    size_t cur_idx = 0;
    while(cur_idx < num_of_elements_to_decode) {
	Frame_t cur_frame;
	cur_frame.read_from_buffer(input_buffer);

	if(cur_frame.num_of_elements > cur_frame.parameters.size()) {

	    for(size_t i = 0; i < cur_frame.num_of_elements; ++i) {
		values[cur_idx + i] =
		    gsl_poly_eval(cur_frame.parameters.data(),
				  cur_frame.parameters.size(),
				  i);
	    }

	} else {

	    std::copy(cur_frame.parameters.begin(),
		      cur_frame.parameters.end(),
		      values.begin() + cur_idx);

	}

	cur_idx += cur_frame.num_of_elements;
    }
}

////////////////////////////////////////////////////////////////////////////////

void
decompress_obt_times(Byte_buffer_t & buffer,
		     double first_obt,
		     size_t num_of_samples,
		     std::vector<uint64_t> & dest)
{
    std::vector<uint32_t> obt_delta_values;
    rle_decompression(buffer, num_of_samples, obt_delta_values);

    dest.resize(obt_delta_values.size() + 1);
    dest[0] = first_obt;

    for(size_t idx = 1; idx < dest.size(); ++idx) {
	dest[idx] = dest[idx - 1] + obt_delta_values[idx - 1];
    }
}

////////////////////////////////////////////////////////////////////////////////

void
decompress_scet_times(Byte_buffer_t & buffer,
		      const Squeezer_file_header_t & file_header,
		      const std::vector<uint64_t> & obt_times,
		      std::vector<double> & dest)
{
    dest.resize(obt_times.size());

    const double slope = 
	(file_header.last_scet_in_ms - file_header.first_scet_in_ms) / 
	(file_header.last_obt - file_header.first_obt);

    for(size_t idx = 0; idx < obt_times.size(); ++idx) {
	double interpolated_scet =
	    file_header.first_scet_in_ms + slope * (obt_times[idx] - file_header.first_obt);
	double scet_correction = buffer.read_float();

	dest[idx] = interpolated_scet + scet_correction;
    }
}

////////////////////////////////////////////////////////////////////////////////

void
decompress_angles(Byte_buffer_t & buffer,
		  size_t num_of_samples,
		  std::vector<double> & dest)
{
    dest.resize(num_of_samples);

    poly_fit_decode(num_of_samples, buffer, dest);

    // Clip angles within [0, 2pi]
    double offset = 0.0;
    for(size_t idx = 0; idx < dest.size(); ++idx) {
	if(dest[idx] + offset < 0.0)
	    offset += 2 * M_PI;
	else if(dest[idx] + offset >= 2 * M_PI)
	    offset -= 2 * M_PI;

	dest[idx] += offset;
    }
}

////////////////////////////////////////////////////////////////////////////////

void
decompress_scientific_data(Byte_buffer_t & buffer,
			   size_t num_of_samples,
			   std::vector<double> & dest)
{
    dest.resize(num_of_samples);
    for(size_t idx = 0; idx < num_of_samples; ++idx)
	dest[idx] = buffer.read_float();
}

////////////////////////////////////////////////////////////////////////////////

void
decompress_quality_flags(Byte_buffer_t & buffer,
			 size_t num_of_samples,
			 std::vector<uint32_t> & dest)
{
    rle_decompression(buffer, num_of_samples, dest);
}

////////////////////////////////////////////////////////////////////////////////

template<typename Process_chunk_data_fn_t>
void
decompress_chunk(size_t chunk_idx,
		 const Squeezer_file_header_t & file_header,
		 const Squeezer_chunk_header_t & chunk_header,
		 std::istream & input_stream,
		 Process_chunk_data_fn_t fn)
{
    if(! chunk_header.is_valid()) {
	throw SqueezerError(boost::format("wrong squeezer chunk %1%")
			    % (chunk_idx + 1));
    }

    Byte_buffer_t chunk_data;
    chunk_data.buffer.resize(chunk_header.number_of_bytes);
    input_stream.read(reinterpret_cast<char *>(chunk_data.buffer.data()),
		      chunk_header.number_of_bytes);
    if(input_stream.bad() || input_stream.fail()) {	
	throw SqueezerError(boost::format("unable to read squeezer chunk %1%")
			    % (chunk_idx + 1));
    }

    fn(chunk_idx, chunk_data, file_header, chunk_header);
}

////////////////////////////////////////////////////////////////////////////////

template<typename Process_chunk_data_fn_t>
void
decompress_file(const std::string & file_name,
		Process_chunk_data_fn_t process_chunk_data_fn,
		Squeezer_file_type_t expected_type)
{
    std::ifstream input_stream(file_name);
    if(input_stream.bad()) {
	throw SqueezerError(boost::format("Unable to open file %1%")
			    % file_name);
    }

    Squeezer_file_header_t file_header(input_stream);
    if(! file_header.is_valid()) {
	throw SqueezerError(boost::format("File %1% does not seem to have been "
					  "created by Squeezer")
			    % file_name);
    }

    if(file_header.get_type() != expected_type) {
	switch(expected_type) {
	case SQZ_DETECTOR_POINTINGS:
	    throw SqueezerError(boost::format("squeezer file %1% does not "
					      "contain pointing information")
				% file_name);
	case SQZ_DIFFERENCED_DATA:
	    throw SqueezerError(boost::format("squeezer file %1% does not "
					      "contain differenced data")
				% file_name);
	default:
	    throw SqueezerError(boost::format("squeezer file %1% has type %2%, but "
					      "type %3% was expected")
				% file_name % file_header.get_type() % expected_type);
	}
    }

    for(size_t idx = 0; idx < file_header.number_of_chunks; ++idx) {
	Squeezer_chunk_header_t chunk_header(input_stream);

	decompress_chunk(idx, file_header, chunk_header, 
			 input_stream, process_chunk_data_fn);
    }
}

////////////////////////////////////////////////////////////////////////////////

struct Process_pointing_chunk {
    PointingData & pnt;
    std::string file_name;

    Process_pointing_chunk(PointingData & a_pnt, std::string a_file_name)
	: pnt(a_pnt), file_name(a_file_name) { }

    void operator()(size_t chunk_idx,
		    Byte_buffer_t & chunk_data,
		    const Squeezer_file_header_t & file_header,
		    const Squeezer_chunk_header_t & chunk_header) 
    {
	switch(chunk_header.chunk_type) {
	case CHUNK_DELTA_OBT: {
	    decompress_obt_times(chunk_data, file_header.first_obt,
				 chunk_header.number_of_samples,
				 pnt.obt_time);
	    break;
	}
	case CHUNK_SCET_ERROR: {
	    if(pnt.obt_time.empty()) {
		throw SqueezerError(boost::format("malformed chunk %1% in file %2%: "
						  "SCET times found before OBT times")
				    % (chunk_idx + 1) % file_name);
	    }

	    decompress_scet_times(chunk_data, file_header, pnt.obt_time, pnt.scet_time);
	    break;
	}
	case CHUNK_THETA: {
	    decompress_angles(chunk_data, chunk_header.number_of_samples, pnt.theta);
	    break;
	}
	case CHUNK_PHI: {
	    decompress_angles(chunk_data, chunk_header.number_of_samples, pnt.phi);
	    break;
	}
	case CHUNK_PSI:	{
	    decompress_angles(chunk_data, chunk_header.number_of_samples, pnt.psi);
	    break;
	}
	default:
	    throw SqueezerError(boost::format("unexpected chunk type %1% while "
					      "decompressing pointing information "
					      "from file %2%")
				% chunk_header.chunk_type % file_name);
	}
    }
};

void
decompress_pointings(const std::string & file_name,
		     PointingData & pointings)
{
    Process_pointing_chunk process_fn(pointings, file_name);
    decompress_file(file_name, process_fn, SQZ_DETECTOR_POINTINGS);
}

////////////////////////////////////////////////////////////////////////////////

struct Process_datadiff_chunk {
    DifferencedData & datadiff;
    std::string file_name;

    Process_datadiff_chunk(DifferencedData & a_datadiff, std::string a_file_name)
	: datadiff(a_datadiff), file_name(a_file_name) { }

    void operator()(size_t chunk_idx,
		    Byte_buffer_t & chunk_data,
		    const Squeezer_file_header_t & file_header,
		    const Squeezer_chunk_header_t & chunk_header) 
    {
	switch(chunk_header.chunk_type) {
	case CHUNK_DELTA_OBT: {
	    decompress_obt_times(chunk_data, file_header.first_obt,
				 chunk_header.number_of_samples,
				 datadiff.obt_time);
	    break;
	}
	case CHUNK_SCET_ERROR: {
	    if(datadiff.obt_time.empty()) {
		throw SqueezerError(boost::format("malformed chunk %1% in file %2%: "
						  "SCET times found before OBT times")
				    % (chunk_idx + 1) % file_name);
	    }

	    decompress_scet_times(chunk_data, file_header,
				  datadiff.obt_time, datadiff.scet_time);
	    break;
	}
	case CHUNK_DIFFERENCED_DATA:
	{
	    decompress_scientific_data(chunk_data, 
				       chunk_header.number_of_samples,
				       datadiff.sky_load);
	    break;
	}
	case CHUNK_QUALITY_FLAGS:
	{
	    decompress_quality_flags(chunk_data, 
				     chunk_header.number_of_samples,
				     datadiff.flags);
	    break;
	}
	default:
	    throw SqueezerError(boost::format("unexpected chunk type %1% while "
					      "decompressing differenced data "
					      "from file %2%")
				% chunk_header.chunk_type % file_name);
	}
    }
};

void
decompress_differenced_data(const std::string & file_name,
			    DifferencedData & datadiff)
{
    Process_datadiff_chunk process_fn(datadiff, file_name);
    decompress_file(file_name, process_fn, SQZ_DIFFERENCED_DATA);
}
