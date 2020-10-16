    /**
 * \section License
 * This file is part of JuPedSim.
 *
 * JuPedSim is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * JuPedSim is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with JuPedSim. If not, see <http://www.gnu.org/licenses/>.
 **/
//
// Created by laemmel on 23.03.16.
//
#pragma once

#include "JPSfire/B_walking_speed/WalkingSpeed.h"
#include "JPSfire/C_toxicity_analysis/ToxicityAnalysis.h"
#include "Macros.h"
#include "direction/DirectionManager.h"
#include "general/Filesystem.h"
#include "math/OperationalModel.h"
#include "pedestrian/AgentsParameters.h"
#include "randomnumbergenerator.h"
#include "routing/RoutingEngine.h"

#include <cstdlib>
#include <memory>
#include <set>
#include <string>


// This class provides a data container for all configuration parameters.
class Configuration
{
public:
    Configuration()
    {
        _walkingSpeed     = nullptr;
        _ToxicityAnalysis = nullptr;
        _routingEngine    = std::shared_ptr<RoutingEngine>(new RoutingEngine());
        _maxOpenMPThreads = 1;
        _seed             = 0;
        _fps              = 8;
        _precision        = 2;
        _linkedCellSize   = 2.2;     // meter
        _model            = nullptr; // std::shared_ptr<OperationalModel>(new OperationalModel());
        _tMax             = 500;     // seconds
        _dT               = 0.01;
        _isPeriodic       = 0; // use only for Tordeux2015 with "trivial" geometries
        // ----------- GCFM repulsive force ------
        _nuPed  = 0.4;
        _nuWall = 0.2;
        // ----------- repulsive force ------
        _aPed  = 1;    // Tordeux2015
        _bPed  = 0.25; // Tordeux2015
        _cPed  = 3;
        _aWall = 1;
        _bWall = 0.7;
        _cWall = 3;
        // ----------- Tordeux2015 model ------
        _dWall = 0.1;
        _dPed  = 0.1;
        // ------- Interpolation GCFM - left side
        _intPWidthPed  = 0.1;
        _intPWidthWall = 0.1;
        // ------- GCFM repulsive force
        _maxFPed  = 3;
        _maxFWall = 3;
        // -------- Interpolation GCFM - right side
        _distEffMaxPed  = 2;
        _distEffMaxWall = 2;
        // ----------------

        _hostname                 = "localhost";
        _trajectoriesFile         = "trajectories.txt";
        _originalTrajectoriesFile = "trajectories.txt";
        _errorLogFile             = "log.txt";
        _projectFile              = "";
        _geometryFile             = "";
        _projectRootDir           = ".";
        _showStatistics           = false;
        _fileFormat               = FileFormat::TXT;
        _agentsParameters         = std::map<int, std::shared_ptr<AgentsParameters>>();
        // ---------- floorfield
        _deltaH              = 0.0625;
        _wall_avoid_distance = 0.4;
        _use_wall_avoidance  = true;

        // ff router quickest
        _recalc_interval = 3;

        // ff router
        _has_specific_goals         = false;
        _has_directional_escalators = false;
        _write_VTK_files            = false;
        _exit_strat                 = 9;
        _write_VTK_files_direction  = false;
        _dirManager                 = nullptr;
        // for random numbers
        _rdGenerator = RandomNumberGenerator();
    }

    std::shared_ptr<WalkingSpeed> GetWalkingSpeed() { return _walkingSpeed; };
    void SetWalkingSpeed(std::shared_ptr<WalkingSpeed> & w) { _walkingSpeed = w; };

    std::shared_ptr<ToxicityAnalysis> GetToxicityAnalysis() { return _ToxicityAnalysis; };
    void SetToxicityAnalysis(std::shared_ptr<ToxicityAnalysis> & t) { _ToxicityAnalysis = t; };

    std::shared_ptr<RoutingEngine> GetRoutingEngine() const { return _routingEngine; };

    // TODO: this is certainly not a config parameter but part of the model, we
    // really should separate data and model [gl march '16]
    void SetRoutingEngine(std::shared_ptr<RoutingEngine> routingEngine)
    {
        _routingEngine = routingEngine;
    };

    int GetMaxOpenMPThreads() const { return _maxOpenMPThreads; };

    void SetMaxOpenMPThreads(int maxOpenMPThreads) { _maxOpenMPThreads = maxOpenMPThreads; };

    unsigned int GetSeed() const { return _seed; };

    void SetSeed(unsigned int seed) { _seed = seed; };

    double GetFps() const { return _fps; };

    void SetFps(double fps) { _fps = fps; };

    unsigned int GetPrecision() const { return _precision; };

    void SetPrecision(unsigned int precision) { _precision = precision; };

    double GetLinkedCellSize() const { return _linkedCellSize; };

    void SetLinkedCellSize(double linkedCellSize) { _linkedCellSize = linkedCellSize; };

    std::shared_ptr<OperationalModel> GetModel() const { return _model; };

    void SetModel(std::shared_ptr<OperationalModel> model) { _model = model; };

    double GetTmax() const { return _tMax; };

    void SetTmax(double tMax) { _tMax = tMax; };

    double Getdt() const { return _dT; };

    void Setdt(double dT) { _dT = dT; };

    int IsPeriodic() const { return _isPeriodic; };

    void SetIsPeriodic(int isPeriodic) { _isPeriodic = isPeriodic; };

    double GetNuPed() const { return _nuPed; };

    void SetNuPed(double nuPed) { _nuPed = nuPed; };

    double GetNuWall() const { return _nuWall; };

    void SetNuWall(double nuWall) { _nuWall = nuWall; };

    double GetaPed() const { return _aPed; };

    void SetaPed(double aPed) { _aPed = aPed; };

    double GetbPed() const { return _bPed; };

    void SetbPed(double bPed) { _bPed = bPed; };

    double GetcPed() const { return _cPed; };

    void SetcPed(double cPed) { _cPed = cPed; };

    double GetaWall() const { return _aWall; };

    void SetaWall(double aWall) { _aWall = aWall; };

    double GetbWall() const { return _bWall; };

    void SetbWall(double bWall) { _bWall = bWall; };

    double GetcWall() const { return _cWall; };

    void SetcWall(double cWall) { _cWall = cWall; };

    double GetDWall() const { return _dWall; };

    void SetDWall(double dWall) { _dWall = dWall; };

    double GetDPed() const { return _dPed; };

    void SetDPed(double dPed) { _dPed = dPed; };

    double GetIntPWidthPed() const { return _intPWidthPed; };

    void SetIntPWidthPed(double intPWidthPed) { _intPWidthPed = intPWidthPed; };

    double GetIntPWidthWall() const { return _intPWidthWall; };

    void SetIntPWidthWall(double intPWidthWall) { _intPWidthWall = intPWidthWall; };

    double GetMaxFPed() const { return _maxFPed; };

    void SetMaxFPed(double maxFPed) { _maxFPed = maxFPed; };

    double GetMaxFWall() const { return _maxFWall; };

    void SetMaxFWall(double maxFWall) { _maxFWall = maxFWall; };

    double GetDistEffMaxPed() const { return _distEffMaxPed; };

    void SetDistEffMaxPed(double distEffMaxPed) { _distEffMaxPed = distEffMaxPed; };

    double GetDistEffMaxWall() const { return _distEffMaxWall; };
    
    double get_xmax() const { return _xmax; };

    void set_xmax(double xmax) { _xmax = xmax; };

    double get_ymin() const { return _ymin; };

    void set_ymin(double ymin) { _ymin = ymin; };

    double get_ymax() const { return _ymax; };

    void set_ymax(double ymax) { _ymax = ymax; };

    double get_xmin() const { return _xmin; };

    void set_xmin(double xmin) { _xmin = xmin; };

    double get_emu() const { return _emu;};

    void set_emu(double emu) { _emu = emu;};

    double get_esigma() const { return _esigma;};

    void set_esigma(double esigma) { _esigma = esigma;};
    
    void SetDistEffMaxWall(double distEffMaxWall) { _distEffMaxWall = distEffMaxWall; };

    double get_deltaH() const { return _deltaH; }

    void set_deltaH(double deltaH) { _deltaH = deltaH; }

    double get_wall_avoid_distance() const { return _wall_avoid_distance; }

    void set_wall_avoid_distance(double wall_avoid_distance)
    {
        _wall_avoid_distance = wall_avoid_distance;
    }

    bool get_use_wall_avoidance() const { return _use_wall_avoidance; }

    void set_use_wall_avoidance(bool use_wall_avoidance)
    {
        _use_wall_avoidance = use_wall_avoidance;
    }

    double get_recalc_interval() const { return _recalc_interval; }

    void set_recalc_interval(double recalc_interval) { _recalc_interval = recalc_interval; }

    bool get_has_specific_goals() const { return _has_specific_goals; }

    void set_has_specific_goals(bool has_specific_goals)
    {
        _has_specific_goals = has_specific_goals;
    }

    bool get_has_directional_escalators() const { return _has_directional_escalators; }
    void set_has_directional_escalators(bool has_directional_esc)
    {
        _has_directional_escalators = has_directional_esc;
    }

    void set_write_VTK_files(bool write_VTK_files) { _write_VTK_files = write_VTK_files; }

    bool get_write_VTK_files() const { return _write_VTK_files; }

    void set_exit_strat(int e_strat) { _exit_strat = e_strat; }

    int get_exit_strat() const { return _exit_strat; }

    void SetDirectionManager(std::shared_ptr<DirectionManager> dir) { _dirManager = dir; }
    std::shared_ptr<DirectionManager> GetDirectionManager() { return _dirManager; }

    const std::string & GetHostname() const { return _hostname; };

    void set_write_VTK_files_direction(bool write_VTK_files_direction)
    {
        _write_VTK_files_direction = write_VTK_files_direction;
    }

    bool get_write_VTK_files_direction() const { return _write_VTK_files_direction; }

    void SetHostname(std::string hostname) { _hostname = hostname; };

    const fs::path & GetTrajectoriesFile() const { return _trajectoriesFile; };

    void SetTrajectoriesFile(const fs::path & trajectoriesFile)
    {
        _trajectoriesFile = trajectoriesFile;
    };

    const fs::path & GetOutputPath() const { return _outputPath; };

    void SetOutputPath(const fs::path & outputDir) { _outputPath = outputDir; };

    void ConfigureOutputPath()
    {
        // Set default name if none was set
        if(_outputPath.empty()) {
            _outputPath = "results";
        }

        // make absolute path
        if(_outputPath.is_relative()) {
            _outputPath = fs::absolute(_outputPath);
        }

        // checks if directory exists, otherwise creates it.
        fs::create_directories(_outputPath);
    }

    const fs::path & GetOriginalTrajectoriesFile() const { return _originalTrajectoriesFile; };

    void SetOriginalTrajectoriesFile(const fs::path & trajectoriesFile)
    {
        _originalTrajectoriesFile = trajectoriesFile;
    };

    const fs::path & GetErrorLogFile() const { return _errorLogFile; };

    void SetErrorLogFile(const fs::path & errorLogFile) { _errorLogFile = errorLogFile; };

    const fs::path & GetProjectFile() const { return _projectFile; };

    void SetProjectFile(const fs::path & projectFile) { _projectFile = projectFile; };

    const fs::path & GetGeometryFile() const { return _geometryFile; };

    void SetGeometryFile(const fs::path & geometryFile) { _geometryFile = geometryFile; };

    const fs::path & GetTransitionFile() const { return _transitionFile; }

    void SetTransitionFile(const fs::path & transitionFile) { _transitionFile = transitionFile; }

    const fs::path & GetGoalFile() const { return _goalFile; }

    void SetGoalFile(const fs::path & goalFile) { _goalFile = goalFile; }

    const fs::path & GetSourceFile() const { return _sourceFile; }

    void SetSourceFile(const fs::path & sourceFile) { _sourceFile = sourceFile; }

    const fs::path & GetTrafficContraintFile() const { return _trafficContraintFile; }

    void SetTrafficContraintFile(const fs::path & trafficContraintFile)
    {
        _trafficContraintFile = trafficContraintFile;
    }

    const fs::path & GetEventFile() const { return _eventFile; }

    void SetEventFile(const fs::path & eventFile) { _eventFile = eventFile; }

    const fs::path & GetScheduleFile() const { return _scheduleFile; }

    void SetScheduleFile(const fs::path & scheduleFile) { _scheduleFile = scheduleFile; }

    const fs::path & GetTrainTypeFile() const { return _trainTypeFile; }
    void SetTrainTypeFile(const fs::path & trainTypeFile) { _trainTypeFile = trainTypeFile; }
    const fs::path & GetTrainTimeTableFile() const { return _trainTimeTableFile; }
    void SetTrainTimeTableFile(const fs::path & trainTimeTableFile)
    {
        _trainTimeTableFile = trainTimeTableFile;
    }

    const fs::path & GetProjectRootDir() const { return _projectRootDir; };

    void SetProjectRootDir(const fs::path & projectRootDir) { _projectRootDir = projectRootDir; };

    bool ShowStatistics() const { return _showStatistics; };

    void SetShowStatistics(bool showStatistics) { _showStatistics = showStatistics; };

    const FileFormat & GetFileFormat() const { return _fileFormat; };

    void SetFileFormat(FileFormat fileFormat) { _fileFormat = fileFormat; };

    const std::map<int, std::shared_ptr<AgentsParameters>> & GetAgentsParameters() const
    {
        return _agentsParameters;
    };

    void AddAgentsParameters(std::shared_ptr<AgentsParameters> agentsParameters, int id)
    {
        _agentsParameters[id] = agentsParameters;
    };

    RandomNumberGenerator * GetRandomNumberGenerator() const { return &_rdGenerator; };

    void AddOptionalOutputOption(OptionalOutput option) { _optionalOutput.insert(option); };

    std::set<OptionalOutput> GetOptionalOutputOptions() { return _optionalOutput; };

private:
    std::shared_ptr<WalkingSpeed> _walkingSpeed;
    std::shared_ptr<ToxicityAnalysis> _ToxicityAnalysis;
    std::shared_ptr<RoutingEngine> _routingEngine;
    int _maxOpenMPThreads;
    unsigned int _seed;
    double _fps;
    unsigned int _precision;
    double _linkedCellSize;
    std::shared_ptr<OperationalModel> _model;
    double _tMax;
    double _dT;
    int _isPeriodic;
    double _nuPed;
    double _nuWall;
    double _aPed;
    double _bPed;
    double _cPed;
    double _aWall;
    double _bWall;
    double _cWall;
    double _dWall;
    double _dPed;
    double _xmin;
    double _ymin;
    double _xmax;
    double _ymax;
    double _emu;
    double _esigma;
    double _intPWidthPed;
    double _intPWidthWall;
    double _maxFPed;
    double _maxFWall;
    double _distEffMaxPed;
    double _distEffMaxWall;
    // floorfield
    double _deltaH;
    double _wall_avoid_distance;
    bool _use_wall_avoidance;

    // ff router quickest
    double _recalc_interval;

    // ff router
    bool _has_specific_goals;
    bool _has_directional_escalators;
    bool _write_VTK_files;
    bool _write_VTK_files_direction;

    int _exit_strat;

    std::shared_ptr<DirectionManager> _dirManager;

    std::string _hostname;
    fs::path _trajectoriesFile;
    fs::path _originalTrajectoriesFile;
    fs::path _errorLogFile;
    fs::path _projectFile;
    fs::path _geometryFile;
    fs::path _transitionFile;
    fs::path _goalFile;
    fs::path _sourceFile;
    fs::path _trafficContraintFile;
    fs::path _eventFile;
    fs::path _scheduleFile;
    fs::path _trainTypeFile;
    fs::path _trainTimeTableFile;

private:
    fs::path _projectRootDir;
    fs::path _outputPath;
    bool _showStatistics;

    mutable RandomNumberGenerator _rdGenerator;

    FileFormat _fileFormat;
    std::map<int, std::shared_ptr<AgentsParameters>> _agentsParameters;

    std::set<OptionalOutput> _optionalOutput;
};
