#!/usr/bin/env cython
#-*- coding:utf-8 -*-

""" Generation tools for NNGT """

from libc.stdint cimport uintptr_t

from libcpp cimport bool

from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector


# ---------------------- #
# Load the c++ functions #
# ---------------------- #

cdef extern from "../libgrowth/elements_types.hpp" namespace "growth":
    ctypedef pair[vector[double], vector[double]] SkelNeurite
    ctypedef vector[vector[double]] SkelSomas


cdef extern from "../libgrowth/config.hpp" namespace "growth":
    ctypedef enum dtype: BOOL, DOUBLE, INT, VEC_LONG, STRING

    ctypedef unordered_map[string, Property] statusMap

    cdef cppclass Property:
        Property()
        Property(const vector[long]& v)
        Property(const string& s)
        dtype data_type
        bool b
        double d
        int i
        vector[long] l
        string s


cdef extern from "../libgrowth/growth_time.hpp" namespace "growth":
    cdef cppclass Time:
        Time()
        Time(float seconds, unsigned char minutes, unsigned char hours,
             unsigned char days)
        Time(Time initial_time, unsigned long steps)
        void set_sec(float seconds)
        void set_min(char minutes)
        void set_hour(char hours)
        void set_day(char days)
        float get_sec() const
        char get_min() const
        char get_hour() const
        char get_day() const
        double get_total_seconds() const


# exception
cdef extern from "../libgrowth/exceptions.hpp" namespace "growth":
    cdef cppclass BadPropertyType:
        BadPropertyType(string&, string&, string&) except +


# kernel functions
ctypedef unordered_map[ string, vector[double] ] mapParams


cdef extern from "../module.hpp" namespace "growth":
    cdef void init_growth( int* argc, char** argv[] ) except +
    cdef void finalize_growth() except +

    cdef size_t create_objects(const string& object_name,
                               const vector[statusMap]& obj_params) except +
    cdef size_t create_neurons(const vector[statusMap]& neuron_params,
                               const vector[statusMap]& axon_params,
                               const vector[statusMap]& dendrites_params
                               ) except +

    cdef void get_environment ( GEOSGeometry* &environment ) except +

    cdef void set_environment ( GEOSGeometry* environment ) except +

    cdef const Time get_current_time() except +

    cdef statusMap get_kernel_status() except +

    cdef size_t get_num_objects() except +

    cdef void get_skeleton(
        SkelNeurite& axon, SkelNeurite& dendrites, SkelNeurite& nodes,
        SkelNeurite& growth_cones, SkelSomas& somas,
        vector[size_t] gids) except +

    cdef void get_swc(string output_file,
        vector[size_t] gids, unsigned int resolution) except +

    cdef statusMap get_status(size_t gid) except +
    cdef statusMap get_neurite_status(size_t gid,
                                      const string& n_type) except +

    cdef vector[size_t] get_neurons() except +

    cdef void get_defaults(const string& object_name,
                           const string& object_type, statusMap &status) except +
    cdef void get_models(vector[string]& models, const string& object_type) except +

    cdef string object_type(size_t gid) except +

    cdef void reset_kernel() except +

    cdef void set_kernel_status(statusMap status_dict, string c_simulation_ID) except +

    cdef string get_simulation_ID() except +

    cdef void set_status(size_t gid, statusMap neuron_status,
                                     statusMap axon_status,
                                     statusMap dendrites_status) except +

    cdef void simulate(const Time& simtime) except +

    cdef void test_random_generator(vector[vector[double]]& values, size_t size) except+


# ---------------------- #
# GEOS-related functions #
# ---------------------- #

cdef extern from "geos_c.h":
    ctypedef struct GEOSGeometry


cdef inline GEOSGeometry *geos_from_shapely(shapely_geom) except *:
    '''
    Get the GEOS geometry pointer from the given shapely geometry.
    '''
    cdef uintptr_t geos_geom = shapely_geom._geom
    return <GEOSGeometry *>geos_geom
