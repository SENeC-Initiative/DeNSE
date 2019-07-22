/*
 * recorders.hpp
 *
 * This file is part of DeNSE.
 *
 * Copyright (C) 2019 SeNEC Initiative
 *
 * DeNSE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * DeNSE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DeNSE. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef RECORDERS_H
#define RECORDERS_H

// C++ includes
#include <array>
#include <fstream>
#include <iostream>
#include <unordered_set>

// includes from libgrowth
#include "config.hpp"
#include "elements_types.hpp"
#include "growth_time.hpp"


namespace growth
{

/*
 * Typedefs for recordings
 */
typedef void (*ptRecorder)(GCPtr, std::ofstream &);

typedef std::unordered_map<std::string, std::vector<double>> mapNameVecDouble;
typedef std::unordered_map<std::string, std::vector<Time>> mapNameVecTime;
typedef std::unordered_map<stype, std::vector<double>> mapNumVecDouble;
typedef std::unordered_map<stype, std::array<Time, 2>> mapNumArrayTime;

typedef std::unordered_map<stype, std::vector<double>> neuronRec;
typedef std::unordered_map<stype, std::vector<Time>> neuronTimeRec;

typedef std::unordered_map<stype, mapNameVecDouble> neuriteRec;
typedef std::unordered_map<stype, mapNameVecTime> neuriteTimeRec;

typedef std::unordered_map<stype,
                           std::unordered_map<std::string, mapNumVecDouble>>
    gcRec;
typedef std::unordered_map<stype,
                           std::unordered_map<std::string, mapNumArrayTime>>
    gcContinuousTimes;
typedef std::unordered_map<
    stype, std::unordered_map<std::string,
                               std::unordered_map<stype, std::vector<Time>>>>
    gcDiscreteTimes;
typedef std::unordered_map<
    stype, std::unordered_map<std::string, std::unordered_set<stype>>>
      setDeadCones;
typedef std::unordered_map<
    stype, std::unordered_map<std::string, std::unordered_map<stype, stype>>>
    gcNumTimes;


/*
 * Base class for all recorders
 *
 * @warning a Recorder records only from neurons on a single OMP thread.
 */
class BaseRecorder
{
  public:
    BaseRecorder();

    virtual void record();
    virtual void record(const Event &ev);

    virtual unsigned int get_event_type() const;
    virtual unsigned int get_level() const;
    virtual bool get_next_recording(std::vector<Property> &ids,
                                    std::vector<double> &values);
    virtual bool get_next_time(std::vector<Property> &ids,
                               std::vector<double> &values,
                               const std::string &time_units);
    void reset_iterations();
    virtual void final_timestep(stype step);

    void neuron_deleted(stype gid);
    virtual void new_neurite(stype neuron, const std::string& neurite) {};
    virtual void gc_died(stype neuron, const std::string& neurite,
                         stype gc_id) {};

    void get_status(statusMap &status) const;
    virtual void set_status(const statusMap &status);

  protected:
    // state parameters
    std::string record_to_;     // record destination (memory or filename)
    std::ofstream record_file_; // access the file
    int interval_;              // interval between successive records
    std::string restrict_to_;   // dendrite name
    // ids and pointers to targets
    std::unordered_map<stype, NeuronPtr> targets_;
    bool v_iterating_;
    bool t_iterating_;
    // model-related parameters
    std::string observable_;
    const char *cstr_obs_;
};


/*
 * Recording at neuron level
 */
class NeuronContinuousRecorder : public BaseRecorder
{
  public:
    NeuronContinuousRecorder();

    virtual void record() override;
    virtual void final_timestep(stype step) override;

    virtual unsigned int get_event_type() const override;
    virtual unsigned int get_level() const override;
    virtual bool get_next_recording(std::vector<Property> &ids,
                                    std::vector<double> &values) override;
    virtual bool get_next_time(std::vector<Property> &ids,
                               std::vector<double> &values,
                               const std::string &time_units) override;
    virtual void set_status(const statusMap &status) override;

  private:
    neuronRec recording_;
    std::array<Time, 2> times_;
    stype num_times_;
    neuronRec::const_iterator neuron_it_;
};


class NeuronDiscreteRecorder : public BaseRecorder
{
  public:
    NeuronDiscreteRecorder();

    virtual void record(const Event &ev) override;

    virtual unsigned int get_event_type() const override;
    virtual unsigned int get_level() const override;
    virtual bool get_next_recording(std::vector<Property> &ids,
                                    std::vector<double> &values) override;
    virtual bool get_next_time(std::vector<Property> &ids,
                               std::vector<double> &values,
                               const std::string &time_units) override;
    virtual void set_status(const statusMap &status) override;

  private:
    neuronRec recording_;
    neuronTimeRec times_;
    neuronRec::const_iterator neuron_it_;
    neuronTimeRec::const_iterator time_it_;
};


class NeuriteContinuousRecorder : public BaseRecorder
{
  public:
    NeuriteContinuousRecorder();

    virtual void record() override;
    virtual void final_timestep(stype step) override;
    virtual void new_neurite(stype neuron,
                             const std::string& neurite) override;

    virtual unsigned int get_event_type() const override;
    virtual unsigned int get_level() const override;
    virtual bool get_next_recording(std::vector<Property> &ids,
                                    std::vector<double> &values) override;
    virtual bool get_next_time(std::vector<Property> &ids,
                               std::vector<double> &values,
                               const std::string &time_units) override;
    virtual void set_status(const statusMap &status) override;

  private:
    neuriteRec recording_;
    std::array<Time, 2> times_;
    stype num_times_;
    neuriteRec::const_iterator neuron_it_;
    mapNameVecDouble::const_iterator neurite_it_;
    mapNameVecDouble::const_iterator neurite_endit_;
};


class NeuriteDiscreteRecorder : public BaseRecorder
{
  public:
    NeuriteDiscreteRecorder();

    virtual void record(const Event &ev) override;
    virtual void new_neurite(stype neuron,
                             const std::string& neurite) override;

    virtual unsigned int get_event_type() const override;
    virtual unsigned int get_level() const override;
    virtual bool get_next_recording(std::vector<Property> &ids,
                                    std::vector<double> &values) override;
    virtual bool get_next_time(std::vector<Property> &ids,
                               std::vector<double> &values,
                               const std::string &time_units) override;
    virtual void set_status(const statusMap &status) override;

  private:
    neuriteRec recording_;
    neuriteTimeRec times_;
    neuriteRec::const_iterator v_neuron_it_;
    mapNameVecDouble::const_iterator v_neurite_it_;
    mapNameVecDouble::const_iterator v_neurite_endit_;
    neuriteTimeRec::const_iterator t_neuron_it_;
    mapNameVecTime::const_iterator t_neurite_it_;
    mapNameVecTime::const_iterator t_neurite_endit_;
};


class GrowthConeContinuousRecorder : public BaseRecorder
{
  public:
    GrowthConeContinuousRecorder();

    virtual void record() override;
    virtual void record(const Event &ev) override;
    virtual void final_timestep(stype step) override;

    virtual void new_neurite(stype neuron,
                             const std::string& neurite) override;
    virtual void gc_died(stype neuron, const std::string& neurite,
                         stype gc_id) override;

    virtual unsigned int get_event_type() const override;
    virtual unsigned int get_level() const override;
    virtual bool get_next_recording(std::vector<Property> &ids,
                                    std::vector<double> &values) override;
    virtual bool get_next_time(std::vector<Property> &ids,
                               std::vector<double> &values,
                               const std::string &time_units) override;
    virtual void set_status(const statusMap &status) override;

  private:
    gcRec recording_;
    gcContinuousTimes times_;
    gcNumTimes num_times_;
    setDeadCones dead_cones_;
    // value iterators
    gcRec::const_iterator v_neuron_it_;
    std::unordered_map<std::string, mapNumVecDouble>::const_iterator
        v_neurite_it_;
    std::unordered_map<std::string, mapNumVecDouble>::const_iterator
        v_neurite_endit_;
    mapNumVecDouble::const_iterator v_gc_pos_;
    mapNumVecDouble::const_iterator v_gc_endpos_;
    // time iterators
    gcContinuousTimes::const_iterator t_neuron_it_;
    std::unordered_map<std::string, mapNumArrayTime>::const_iterator
        t_neurite_it_;
    std::unordered_map<std::string, mapNumArrayTime>::const_iterator
        t_neurite_endit_;
    mapNumArrayTime::const_iterator t_gc_pos_;
    mapNumArrayTime::const_iterator t_gc_endpos_;
};


class GrowthConeDiscreteRecorder : public BaseRecorder
{
  public:
    GrowthConeDiscreteRecorder();

    virtual void record(const Event &ev) override;

    virtual unsigned int get_event_type() const override;
    virtual unsigned int get_level() const override;
    virtual bool get_next_recording(std::vector<Property> &ids,
                                    std::vector<double> &values) override;
    virtual bool get_next_time(std::vector<Property> &ids,
                               std::vector<double> &values,
                               const std::string &time_units) override;
    virtual void set_status(const statusMap &status) override;

  private:
    gcRec recording_;
    gcDiscreteTimes times_;
    // value iterators
    gcRec::const_iterator v_neuron_it_;
    std::unordered_map<std::string, mapNumVecDouble>::const_iterator
        v_neurite_it_;
    std::unordered_map<std::string, mapNumVecDouble>::const_iterator
        v_neurite_endit_;
    mapNumVecDouble::const_iterator v_gc_pos_;
    mapNumVecDouble::const_iterator v_gc_endpos_;
    // time iterators
    gcDiscreteTimes::const_iterator t_neuron_it_;
    std::unordered_map<std::string,
                       std::unordered_map<stype, std::vector<Time>>>::
        const_iterator t_neurite_it_;
    std::unordered_map<std::string,
                       std::unordered_map<stype, std::vector<Time>>>::
        const_iterator t_neurite_endit_;
    std::unordered_map<stype, std::vector<Time>>::const_iterator t_gc_pos_;
    std::unordered_map<stype, std::vector<Time>>::const_iterator t_gc_endpos_;
};


void _record_length(GCPtr gc, std::ofstream &record_file_);
void _record_CR(GCPtr gc, std::ofstream &record_file_);
void _record_branch_topology(std::string branch_type, stype centrifugal_order,
                             stype branch_size, std::ofstream &record_file_);
void _record_branch_geometry(std::string branch_type, double angle,
                             double ratio, std::ofstream &record_file_);

} // namespace growth

#endif
