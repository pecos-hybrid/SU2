 #pragma once
#include "SzlFileLoader.h"
#include "fileio.h"
namespace tecplot { namespace tecioszl { struct Text : public tecplot::___3933::Text { public: static Text invalidText() { return Text( 0.0, 0.0, 0.0, ___364, ___4455, 0.1, 20.0, ___4075, ___4050, 0.0, ___498, ___364, 14.0, 1.0, "", ___663, ___3446, ___4271, "", "Helvetica", ___1305, ___1305, 1, ___1305 ); } Text( double ___4574, double ___4591, double ___4713, ___516 ___4059, ___516 ___4061, double ___4071, double ___4073, TextBox_e ___4078, TextAnchor_e ___4043, double ___4056, Clipping_e ___4079, ___516 ___4080, double ___4103, double ___4107, std::string ___4109, CoordSys_e ___4115, Scope_e ___4119, Units_e ___4124, std::string ___4126, std::string ___4129, ___372 ___4132, ___372 ___4134, ___1172 ___4138, ___372 ___4105) : tecplot::___3933::Text(___4574, ___4591, ___4713, ___4059, ___4061, ___4071, ___4073, ___4078, ___4043, ___4056, ___4079, ___4080, ___4103, ___4107, ___4109, ___4115, ___4119, ___4124, ___4126, ___4129, ___4132, ___4134, ___4138, ___4105) {} bool ___2067() { return VALID_ENUM(___2639, CoordSys_e); } void writeToFile(std::ofstream& outputFile, bool ___4480) const { writeScalar(outputFile, ___2628, ___4480); writeScalar(outputFile, ___2629, ___4480); writeScalar(outputFile, ___2630, ___4480); writeScalar(outputFile, ___2631, ___4480); writeScalar(outputFile, ___2632, ___4480); writeScalar(outputFile, ___2625, ___4480); writeScalar(outputFile, ___2626[0], ___4480); writeScalar(outputFile, ___2626[1], ___4480); writeScalar(outputFile, ___2626[2], ___4480); writeScalar(outputFile, ___2627, ___4480); writeScalar(outputFile, ___2633, ___4480); writeScalar(outputFile, ___2634, ___4480); writeScalar(outputFile, ___2635, ___4480); writeScalar(outputFile, ___2637, ___4480); ___4544(outputFile, ___2638, ___4480); writeScalar(outputFile, ___2639, ___4480); writeScalar(outputFile, ___2641, ___4480); writeScalar(outputFile, ___2642, ___4480); ___4544(outputFile, ___2643, ___4480); ___4544(outputFile, ___2644, ___4480); writeScalar(outputFile, ___2645, ___4480); writeScalar(outputFile, ___2646, ___4480); writeScalar(outputFile, ___2647, ___4480); writeScalar(outputFile, ___2636, ___4480); } Text(std::ifstream& inputFile, bool readASCII) { readScalar(inputFile, ___2628, readASCII); readScalar(inputFile, ___2629, readASCII); readScalar(inputFile, ___2630, readASCII); readScalar(inputFile, ___2631, readASCII); readScalar(inputFile, (unsigned int&)___2632, readASCII); readScalar(inputFile, (unsigned int&)___2625, readASCII); readScalar(inputFile, ___2626[0], readASCII); readScalar(inputFile, ___2626[1], readASCII); readScalar(inputFile, ___2626[2], readASCII); readScalar(inputFile, ___2627, readASCII); readScalar(inputFile, (unsigned int&)___2633, readASCII); readScalar(inputFile, ___2634, readASCII); readScalar(inputFile, ___2635, readASCII); readScalar(inputFile, ___2637, readASCII); readString(inputFile, ___2638, readASCII); readScalar(inputFile, (unsigned int&)___2639, readASCII); readScalar(inputFile, (unsigned int&)___2641, readASCII); readScalar(inputFile, (unsigned int&)___2642, readASCII); readString(inputFile, ___2643, readASCII); readString(inputFile, ___2644, readASCII); readScalar(inputFile, ___2645, readASCII); readScalar(inputFile, ___2646, readASCII); readScalar(inputFile, ___2647, readASCII); readScalar(inputFile, ___2636, readASCII); } }; }}
