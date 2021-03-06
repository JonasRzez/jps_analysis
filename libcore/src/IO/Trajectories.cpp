#include "Trajectories.h"

#include "geometry/SubRoom.h"
#include "pedestrian/Pedestrian.h"

#include <Logger.h>
#include <tinyxml.h>

static fs::path getTrainTimeTableFileName(const fs::path & projectFile)
{
    fs::path ret{};

    TiXmlDocument doc(projectFile.string());
    if(!doc.LoadFile()) {
        LOG_ERROR("{}", doc.ErrorDesc());
        LOG_ERROR("GetTrainTimeTable could not parse the project file");
        return ret;
    }
    TiXmlNode * xMainNode = doc.RootElement();
    if(xMainNode->FirstChild("train_constraints")) {
        TiXmlNode * xFileNode =
            xMainNode->FirstChild("train_constraints")->FirstChild("train_time_table");

        if(xFileNode) {
            ret = xFileNode->FirstChild()->ValueStr();
        }
        LOG_INFO("train_time_table <{}>", ret.string());
    } else {
        LOG_INFO("No ttt file found");
        return ret;
    }
    return ret;
}

static fs::path getTrainTypeFileName(const fs::path & projectFile)
{
    fs::path ret{};

    TiXmlDocument doc(projectFile.string());
    if(!doc.LoadFile()) {
        LOG_ERROR("{}", doc.ErrorDesc());
        LOG_ERROR("GetTrainType could not parse the project file");
        return ret;
    }
    TiXmlNode * xMainNode = doc.RootElement();
    if(xMainNode->FirstChild("train_constraints")) {
        auto xFileNode = xMainNode->FirstChild("train_constraints")->FirstChild("train_types");
        if(xFileNode)
            ret = xFileNode->FirstChild()->ValueStr();
        LOG_INFO("train_types <{}>", ret.string());
    } else {
        LOG_INFO("No train types file found");
        return ret;
    }
    return ret;
}

/**
 * TXT format implementation
 */
TrajectoriesTXT::TrajectoriesTXT(unsigned int precision) : Trajectories{precision}
{
    // Add header, info and output for speed
    _optionalOutputHeader[OptionalOutput::speed] = "V\t";
    _optionalOutputInfo[OptionalOutput::speed]   = "#V: speed of the pedestrian (in m/s)\n";
    _optionalOutput[OptionalOutput::speed]       = [](const Pedestrian * ped) {
        return fmt::format(FMT_STRING("{:.2f}\t"), ped->GetV().Norm());
    };

    // Add header, info and output for velocity
    _optionalOutputHeader[OptionalOutput::velocity] = "Vx\tVy\t";
    _optionalOutputInfo[OptionalOutput::velocity] =
        "#Vx: x component of the pedestrian's velocity\n"
        "#Vy: y component of the pedestrian's velocity\n";
    _optionalOutput[OptionalOutput::velocity] = [](const Pedestrian * ped) {
        return fmt::format(FMT_STRING("{:.2f}\t{:.2f}\t"), ped->GetV()._x, ped->GetV()._y);
    };

    // Add header, info and output for final_goal
    _optionalOutputHeader[OptionalOutput::final_goal] = "FG\t";
    _optionalOutputInfo[OptionalOutput::final_goal]   = "#FG: id of final goal\n";
    _optionalOutput[OptionalOutput::final_goal]       = [](const Pedestrian * ped) {
        return fmt::format(FMT_STRING("{}\t"), ped->GetFinalDestination());
    };

    // Add header, info and output for intermediate_goal
    _optionalOutputHeader[OptionalOutput::intermediate_goal] = "CG\t";
    _optionalOutputInfo[OptionalOutput::intermediate_goal]   = "#CG: id of current goal\n";
    _optionalOutput[OptionalOutput::intermediate_goal]       = [](const Pedestrian * ped) {
        return fmt::format(FMT_STRING("{}\t"), ped->GetExitIndex());
    };

    // Add header, info and output for desired direction
    _optionalOutputHeader[OptionalOutput::desired_direction] = "Dx\tDy\t";
    _optionalOutputInfo[OptionalOutput::desired_direction] =
        "#Dx: x component of the pedestrian's desired direction\n"
        "#Dy: y component of the pedestrian's desired direction\n";
    _optionalOutput[OptionalOutput::desired_direction] = [](const Pedestrian * ped) {
        return fmt::format(
            FMT_STRING("{:.2f}\t{:.2f}\t"), ped->GetLastE0()._x, ped->GetLastE0()._y);
    };

    // Add header, info and output for spotlight
    _optionalOutputHeader[OptionalOutput::spotlight] = "SPOT\t";
    _optionalOutputInfo[OptionalOutput::spotlight]   = "#SPOT: ped is highlighted\n";
    _optionalOutput[OptionalOutput::spotlight]       = [](Pedestrian * ped) {
        return fmt::format(FMT_STRING("{}\t"), (int) ped->GetSpotlight());
    };

    // Add header, info and output for router
    _optionalOutputHeader[OptionalOutput::router] = "ROUTER\t";
    _optionalOutputInfo[OptionalOutput::router] =
        "#ROUTER: routing strategy used during simulation\n";
    _optionalOutput[OptionalOutput::router] = [](const Pedestrian * ped) {
        return fmt::format(FMT_STRING("{}\t"), ped->GetRoutingStrategy());
    };

    // Add header, info and output for group
    _optionalOutputHeader[OptionalOutput::group] = "GROUP\t";
    _optionalOutputInfo[OptionalOutput::group]   = "#GROUP: group of the pedestrian\n";
    _optionalOutput[OptionalOutput::group]       = [](const Pedestrian * ped) {
        return fmt::format(FMT_STRING("{}\t"), ped->GetGroup());
    };
}

void TrajectoriesTXT::WriteHeader(long nPeds, double fps, Building * building, int, int count)
{
    const fs::path & projRoot(building->GetProjectRootDir());
    const fs::path tmpGeo = projRoot / building->GetGeometryFilename();

    std::string header = fmt::format("#description: jpscore ({:s})\n", JPSCORE_VERSION);
    header.append(fmt::format("#agents: {:d}\n", nPeds));
    header.append(fmt::format("#count: {:d}\n", count));
    header.append(fmt::format("#framerate: {:0.2f}\n", fps));
    header.append(fmt::format(
        "#geometry: {:s}\n", building->GetConfig()->GetGeometryFile().filename().string()));

    // if used: add source file name
    if(!building->GetConfig()->GetSourceFile().empty()) {
        header.append(fmt::format(
            "#sources: {:s}\n", building->GetConfig()->GetSourceFile().filename().string()));
    }

    // if used: add goal file name
    if(!building->GetConfig()->GetGoalFile().empty()) {
        header.append(fmt::format(
            "#goals: {:s}\n", building->GetConfig()->GetGoalFile().filename().string()));
    }

    // if used: add event file name
    if(!building->GetConfig()->GetEventFile().empty()) {
        header.append(fmt::format(
            "#events: {:s}\n", building->GetConfig()->GetEventFile().filename().string()));
    }

    // if used: add trainTimeTable file name
    if(const fs::path trainTimeTableFileName =
           getTrainTimeTableFileName(building->GetProjectFilename());
       !trainTimeTableFileName.empty()) {
        const fs::path tmpTTT = projRoot / trainTimeTableFileName;
        header.append(fmt::format("#trainTimeTable: {:s}\n", tmpTTT.string()));
    }

    // if used: add trainType file name
    if(const fs::path trainTypeFileName = getTrainTypeFileName(building->GetProjectFilename());
       !trainTypeFileName.empty()) {
        const fs::path tmpTT = projRoot / trainTypeFileName;
        header.append(fmt::format("#trainType: {:s}\n", tmpTT.string()));
    }

    header.append("#ID: the agent ID\n");
    header.append("#FR: the current frame\n");
    header.append("#X,Y,Z: the agents coordinates (in metres)\n");
    header.append("#A, B: semi-axes of the ellipse\n");
    header.append("#ANGLE: orientation of the ellipse\n");
    header.append("#COLOR: color of the ellipse\n");

    // Add info for optional output options
    for(const auto & option : _optionalOutputOptions) {
        header.append(_optionalOutputInfo[option]);
    }

    header.append("\n#ID\tFR\tX\tY\tANGLE\tCOLOR\tspeed_nn\tANGLE_int_nn\tIntID\tAngleDirNn\tr_a");
    // Add header for optional output options
    for(const auto & option : _optionalOutputOptions) {
        header.append(_optionalOutputHeader[option]);
    }
    Write(header);
}

void TrajectoriesTXT::WriteGeometry(Building *) {}

void TrajectoriesTXT::WriteFrame(int frameNr, Building * building)
{
    const std::vector<Pedestrian *> & allPeds = building->GetAllPedestrians();
    for(auto ped : allPeds) {
        double x       = ped->GetPos()._x;
        double y       = ped->GetPos()._y;
        //double z       = ped->GetElevation();
        int color      = ped->GetColor();
        double b       = ped->GetEllipse().GetBmax();
        //double b       = ped->GetSmallerAxis();
        double phi     = atan2(ped->GetEllipse().GetSinPhi(), ped->GetEllipse().GetCosPhi());
        double RAD2DEG = 180.0 / M_PI;
        double dir_phi = atan2(ped->GetDirNn()._y,ped->GetDirNn()._x);
    
        double speed_nn = ped->GetSpeedNn();
        double angle_nn_int = ped->GetAngleNn();
        int intID = ped->GetIntID();
        int intIDN = ped->GetIntIDN();
        unsigned int precision = GetPrecision();
        std::string frame      = fmt::format(
            "{:d}\t{:d}\t{:0.{}f}\t{:0.{}f}\t{:0.2f}\t{:d}\t{:0.2f}\t{:0.10f}\t{:d}\t{:0.10f}\t{:0.3f}",
            ped->GetID(),
            frameNr,
            x,
            precision,
            y,
            precision,
            phi * RAD2DEG,
            color,
            speed_nn,
            angle_nn_int,
            intID,
            dir_phi,
            b);
        for(const auto & option : _optionalOutputOptions) {
            frame.append(_optionalOutput[option](ped));
        }

        Write(frame);
    }
}
