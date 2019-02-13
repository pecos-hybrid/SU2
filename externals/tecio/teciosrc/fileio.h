 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include "ThirdPartyHeadersEnd.h"
#include "AsciiOutputInfo.h"
#include "FileReaderInterface.h"
#include "FileWriterInterface.h"
namespace tecplot { namespace tecioszl { class FileIOException : public std::runtime_error { public: FileIOException(char const* what) : std::runtime_error(what) {} virtual ~FileIOException() throw() {} }; template <typename T> inline void writeScalar(___3933::FileWriterInterface& outputFile, T ___4298, bool ___4480) { if (___4480) { char ___416[100]; tecplot::___3933::encodeAsciiValue<T, false, 0>(___416, 100, ___4298); size_t len = strlen(___416); ___416[len] = '\r'; ___416[len + 1] = '\n'; if (outputFile.fwrite(___416, 1, len + 2) != len + 2) throw FileIOException("Error writing scalar to file"); } else { if (outputFile.fwrite(&___4298, sizeof(___4298), 1) != 1) throw FileIOException("Error writing scalar to file"); } } template <typename T> inline uint64_t scalarSizeInFile(T  , bool ___4480) { if (___4480) return ___3933::___199<T, false>::size + 2; else return sizeof(T); } template <typename T> inline void readScalar(___3933::___1399& inputFile, T& ___4298, bool readASCII) { if (readASCII) { char ___416[100]; inputFile.fgets(___416, 100); std::istringstream inputStream(___416); inputStream >> ___4298; if (inputStream.fail()) throw FileIOException("Error reading scalar from file"); } else { if (inputFile.fread((char*)&___4298, sizeof(___4298), 1) != 1) throw FileIOException("Error reading scalar from file"); } } inline void ___4544(___3933::FileWriterInterface& outputFile, std::string const& str, bool ___4480) { uint64_t length = str.size(); writeScalar(outputFile, length, ___4480); if (outputFile.fwrite(str.c_str(), 1, str.size()) != str.size()) throw FileIOException("Error writing string length to file"); if (___4480 && outputFile.fwrite("\r\n", 1, 2) != 2) throw FileIOException("Error writing string to file"); } inline uint64_t stringSizeInFile(std::string const& str, bool ___4480) { uint64_t length = str.size(); uint64_t sizeInFile = scalarSizeInFile(length, ___4480) + str.size(); if (___4480) sizeInFile += 2; return sizeInFile; } inline void readString(___3933::___1399& inputFile, std::string& str, bool readASCII) { uint64_t length; readScalar(inputFile, length, readASCII); str.resize((size_t)length); if (readASCII) { char ___416[100]; inputFile.fgets(___416, 100); std::istringstream inputStream(___416); inputStream >> str; if (inputStream.fail()) throw FileIOException("Error reading string length from file"); } else { if (inputFile.fread(&str[0], 1, (size_t)length) != (size_t)length) throw FileIOException("Error reading string from file"); } } template <typename T> inline void writeVector(___3933::FileWriterInterface& outputFile, std::vector<T> const& vec, bool ___4480) { uint64_t length = vec.size(); writeScalar(outputFile, length, ___4480); if (vec.size() > 0) { if (___4480) { BOOST_FOREACH(T ___4298, vec) { writeScalar(outputFile, ___4298, ___4480); } } else { if (outputFile.fwrite(&(vec[0]), sizeof(T), vec.size()) != vec.size()) throw FileIOException("Error writing vector to file"); } } } template <typename T> inline uint64_t vectorSizeInFile(std::vector<T> const& vec, bool ___4480) { uint64_t length = vec.size(); uint64_t sizeInFile = scalarSizeInFile(length, ___4480); if (vec.size() > 0) sizeInFile += length * scalarSizeInFile(vec[0], ___4480); return sizeInFile; } template <typename T> inline void readVector(___3933::___1399& inputFile, std::vector<T>& vec, bool readASCII) { uint64_t length; readScalar(inputFile, length, readASCII); vec.resize((size_t)length); if (length > 0) { if (readASCII) { for(uint64_t i = 0; i < length; ++i) readScalar(inputFile, vec[i], readASCII); } else { if (inputFile.fread((&vec[0]), sizeof(T), (size_t)length) != (size_t)length) throw FileIOException("Error reading vector from file"); } } } template <> inline void writeVector<std::string>(___3933::FileWriterInterface& outputFile, std::vector<std::string> const& vec, bool ___4480) { uint64_t length = vec.size(); writeScalar(outputFile, length, ___4480); for(size_t i = 0; i < vec.size(); ++i) ___4544(outputFile, vec[i], ___4480); } template <> inline uint64_t vectorSizeInFile<std::string>(std::vector<std::string> const& vec, bool ___4480) { uint64_t length = vec.size(); uint64_t sizeInFile = scalarSizeInFile(length, ___4480); for(size_t i = 0; i < vec.size(); ++i) sizeInFile += stringSizeInFile(vec[i], ___4480); return sizeInFile; } template <> inline void readVector<std::string>(___3933::___1399& inputFile, std::vector<std::string>& vec, bool readASCII) { uint64_t length; readScalar(inputFile, length, readASCII); vec.resize((size_t)length); for(size_t i = 0; i < vec.size(); ++i) readString(inputFile, vec[i], readASCII); } template <typename T> inline void writeVectorOfObjects(___3933::FileWriterInterface& outputFile, std::vector<T> const& vec, bool ___4480) { uint64_t length = vec.size(); writeScalar(outputFile, length, ___4480); for(size_t i = 0; i < vec.size(); ++i) vec[i].writeToFile(outputFile, ___4480); } template <typename T> inline uint64_t vectorOfObjectsSizeInFile(std::vector<T> const& vec, bool ___4480)
{ uint64_t length = vec.size(); uint64_t sizeInFile = scalarSizeInFile(length, ___4480); for(size_t i = 0; i < vec.size(); ++i) sizeInFile += vec[i].sizeInFile(___4480); return sizeInFile; } template <typename T> inline void readVectorOfObjects(___3933::___1399& inputFile, std::vector<T>& vec, bool readASCII) { uint64_t length; readScalar(inputFile, length, readASCII); vec.reserve((size_t)length); for(uint64_t i = 0; i < length; ++i) vec.push_back(T(inputFile, readASCII)); } template <typename T> inline void writeVectorOfPtrs(___3933::FileWriterInterface& outputFile, std::vector<boost::shared_ptr<T> > const& vec, bool ___4480) { uint64_t length = vec.size(); writeScalar(outputFile, length, ___4480); for(size_t i = 0; i < vec.size(); ++i) vec[i]->writeToFile(outputFile, ___4480); } template <typename T> inline uint64_t vectorOfPtrsSizeInFile(std::vector<boost::shared_ptr<T> > const& vec, bool ___4480) { uint64_t length = vec.size(); uint64_t sizeInFile = scalarSizeInFile(length, ___4480); for(size_t i = 0; i < vec.size(); ++i) sizeInFile += vec[i]->sizeInFile(___4480); return sizeInFile; } template <typename T> inline void readVectorOfPtrs(___3933::___1399& inputFile, std::vector<boost::shared_ptr<T> >& vec, bool readASCII) { uint64_t length; readScalar(inputFile, length, readASCII); vec.resize(0); vec.reserve((size_t)length); for(uint64_t i = 0; i < length; ++i) vec.push_back(T::makePtr(inputFile, readASCII)); } template <typename T, typename U> inline void writeMapOfScalarsToPtrs(___3933::FileWriterInterface& outputFile, std::map<T, boost::shared_ptr<U> > const& m, bool ___4480) { uint64_t length = m.size(); writeScalar(outputFile, length, ___4480); typedef std::pair<T, boost::shared_ptr<U> > ValuePair; BOOST_FOREACH (ValuePair const& valuePair, m) { writeScalar(outputFile, valuePair.first, ___4480); valuePair.second->writeToFile(outputFile, ___4480); } } template <typename T, typename U> inline uint64_t mapOfScalarsToPtrsSizeInFile(std::map<T, boost::shared_ptr<U> > const& m, bool ___4480) { uint64_t length = m.size(); uint64_t sizeInFile = scalarSizeInFile(length, ___4480); typedef std::pair<T, boost::shared_ptr<U> > ValuePair; BOOST_FOREACH (ValuePair const& valuePair, m) { sizeInFile += scalarSizeInFile(valuePair.first, ___4480); sizeInFile += valuePair.second->sizeInFile(___4480); } return sizeInFile; } template <typename T, typename U> inline void readMapOfScalarsToPtrs(___3933::___1399& inputFile, std::map<T, boost::shared_ptr<U> >& m, bool readASCII) { uint64_t length; readScalar(inputFile, length, readASCII); for(uint64_t i = 0; i < length; ++i) { T key; readScalar(inputFile, key, readASCII); m[key] = U::makePtr(inputFile, readASCII); } } template <typename T, typename U> inline void writeMapOfPairsToObjects(___3933::FileWriterInterface& outputFile, std::map<T, U> const& m, bool ___4480) { uint64_t length = m.size(); writeScalar(outputFile, length, ___4480); typedef std::pair<T, U> ValuePair; BOOST_FOREACH (ValuePair const& valuePair, m) { writeScalar(outputFile, valuePair.first.first, ___4480); writeScalar(outputFile, valuePair.first.second, ___4480); valuePair.second.writeToFile(outputFile, ___4480); } } template <typename T, typename U> inline uint64_t mapOfPairsToObjectsSizeInFile(std::map<T, U> const& m, bool ___4480) { uint64_t length = m.size(); uint64_t sizeInFile = scalarSizeInFile(length, ___4480); typedef std::pair<T, U> ValuePair; BOOST_FOREACH (ValuePair const& valuePair, m) { sizeInFile += scalarSizeInFile(valuePair.first.first, ___4480); sizeInFile += scalarSizeInFile(valuePair.first.second, ___4480); sizeInFile += valuePair.second.sizeInFile(___4480); } return sizeInFile; } template <typename T, typename U> inline void readMapOfPairsToObjects(___3933::___1399& inputFile, std::map<T, U>& m, bool readASCII) { uint64_t length; readScalar(inputFile, length, readASCII); for(uint64_t i = 0; i < length; ++i) { T key; readScalar(inputFile, key.first, readASCII); readScalar(inputFile, key.second, readASCII); m[key] = U(inputFile, readASCII); } } }}
