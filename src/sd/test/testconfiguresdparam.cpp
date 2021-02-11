/*
 * testconfiguresdparam.cpp
 *
 *  Created on: 2019/2/15
 *      Author: hbfrank
 */


#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <fstream>
#include <boost/algorithm/string.hpp>

struct QA
{
    QA() {;}
    QA(std::string qidx, std::string q, std::string a)
            : QIndex(qidx), QString(q), AString(a) {;}
    std::string QIndex;
    std::string QString;
    std::string AString;
};

std::pair<std::string, std::string> parseQ(const std::string qline)
{
    std::pair<std::string, std::string> q("", "");
    std::string trqline = boost::trim_copy(qline);
    if (trqline.size() < 5) {
        return q;
    }
    if (trqline.substr(0, 1) != "[") {
        return q;
    }
    int idxrlim = trqline.find(']', 3);
    if (idxrlim == std::string::npos) {
        return q;
    }
    q.first = boost::trim_copy(trqline.substr(1, idxrlim-1)); // (idxrlim-1)-1+1
    q.second = boost::trim_copy(trqline.substr(idxrlim));
    return q;
}

std::string parseA(const std::string aline)
{
    std::string answer = "";
    std::string traline = boost::trim_copy(aline);
    if (traline.size() < 2) {
        return answer;
    }
    if (traline[0] != '{') {
        return answer;
    }
    int idxrlim = traline.find('}', 1);
    if (idxrlim == std::string::npos) {
        return answer;
    }
    answer = boost::trim_copy(traline.substr(1, idxrlim-1)); // (idxrlim-1)-1+1
    return answer;
}

std::vector<std::string> readLinesFromFile(const std::string & inputFile,
                              const std::vector<std::string> &
                                  commentChars = std::vector<std::string>(),
                              const bool skipEmptyLines = true,
                              const bool skipBlankLines = true)
{
    std::vector<std::string> lines;
    std::ifstream ifs(inputFile);
    std::string line;
    std::string ltrline; // left-trimmed line
    bool isCommented = false;
    while (std::getline(ifs, line)) {
        if (skipEmptyLines && line.empty()) continue;
        if (skipBlankLines && boost::trim_copy(line).empty()) continue;
        isCommented = false;
        if (!commentChars.empty()) {
            ltrline = boost::trim_left_copy(line);
            for (const std::string & cc : commentChars) {
                if (ltrline.substr(0, cc.size()) == cc) {
                    isCommented = true;
                    break;
                }
            }
        }
        if (isCommented) continue;
        lines.push_back(std::string(line));
    }
    ifs.close();
    return lines;
}

std::map<std::string, QA> parseQALines(const std::string & inputFile)
{
    std::vector<std::string> commentchars { "#" };
    std::vector<std::string> inputlines = readLinesFromFile(inputFile, commentchars);
    int numLines = inputlines.size();
    if (numLines == 0) {
        std::cerr << "Cease configuring! Either the configuration file "
                "does not exist, or it is empty." << std::endl;
        exit(1);
    }
    if (numLines % 2 != 0) {
        std::cerr << "Cease configuring! Either the configuration file "
                "is broken, or it has odd parameter line." << std::endl;
        exit(1);
    }
    std::map<std::string, QA> qas;
    for (int i = 0; i < numLines-1; i++) {
        std::string qline = inputlines[i];
        std::pair<std::string, std::string> q = parseQ(qline);
        std::string aline = inputlines[i+1];
        std::string a = parseA(aline);
        QA qa(q.first, q.second, a);
        if (q.first.empty()) {
            continue;
        }
        qas.insert(std::pair<std::string, QA>(qa.QIndex, qa));
    }
    return qas;
}

static std::map<std::string, std::string> QAIndices = {
        // SDInputOutput
        {"StartPDBFile", "Q1.1"},
        {"OutputPDBFile", "Q1.2"},
        {"RMSDRefPDBFile", "Q1.4"},
        {"SDTotalSteps", "Q1.5"},
        {"SDSeed", "Q1.6"},
        // GenChain
        {"SCTypeFrom", "Q2.1"},
        {"RefPDBFile", "Q2.1.A"},
        {"SeqString", "Q2.1.B"},
        {"SCCoordFrom", "Q2.2"},
        // SD
        {"ZZZTrajOutput", "Q1.3"},
        {"Temperatures", "Q3.1"},
        {"AnnealingScheme", "Q3.1"},
        {"FrictionCoeff", "Q3.2"},
        // ForceField
        {"LocalHBWeight", "Q4.1"},
        {"ZZZSCEnergyTermsWeight", "Q4.2"},
        {"RgRestrain", "Q4.3"},
        {"ZZZRgRestrainRegionFrom", "Q4.3.1"},
        {"CoreRegion", "Q4.3.1.1"},
        {"TotalRMSDRestrain", "Q4.4"},
        {"ZZZRMSDRestrainValue", "Q4.4.1"},
        {"RMSDRestrainedRegion", "Q4.4.2"},
        {"HelixRestrains", "Q4.5"},
        {"ExtendRestrains", "Q4.6"},
        {"LoopRegionWeights", "Q4.7"}
};

std::string findRawAnswer(const std::string & optname,
                          const std::map<std::string, QA> & qas)
{
    std::string ans;
    std::string qindex;
    auto qaiter = qas.end();
    qindex = QAIndices.find(optname)->second;
    qaiter = qas.find(qindex);
    if (qaiter != qas.end()) {
        ans = qaiter->second.AString;
    }
    return ans;
}

bool setupPlainOption(const std::string & optname,
           const std::map<std::string, QA> & qas,
           std::map<std::string, std::string> & options)
{
    std::string rawAnswer = findRawAnswer(optname, qas);
    options[optname] = rawAnswer;
    return true;
}

class ParamOptions {
public:
    virtual void setup(std::map<std::string, QA>) = 0;
    bool isSufficient() {
        for (auto & o : requiredOptions_) {
            if (options_.find(o) == options_.end()) {
                return false;
            }
        }
        return true;
    }
    virtual std::vector<std::string> toLines() = 0;
    virtual ~ParamOptions() {;}
protected:
    std::vector<std::string> toLines(std::string paramName) {
        std::vector<std::string> lines;
        lines.push_back("START " + paramName);
        for (auto & o : options_) {
            lines.push_back(o.first + "\t=\t" + o.second);
        }
        lines.push_back("END " + paramName);
        return lines;
    }
protected:
    std::map<std::string, std::string> options_;
    std::vector<std::string> requiredOptions_;
};

class SDIOParamOptions : public ParamOptions {
public:
    SDIOParamOptions() {
        options_ =  {
            {"EchoStartPDBFile", ""},
            {"OutputCoordFile", ""},
            {"SDSeed", "13579"}
        };
        requiredOptions_ = {
            "SDTotalSteps",
            "OutputPDBFile"
        };
    }

    void setup(std::map<std::string, QA> qas) {
        // OutputPDBFile
        ::setupPlainOption("OutputPDBFile", qas, options_);

        // OutputCoordFile
        ::setupPlainOption("OutputCoordFile", qas, options_);

        // SDSeed
        ::setupPlainOption("SDSeed", qas, options_);

        // SDTotalSteps
        ::setupPlainOption("SDTotalSteps", qas, options_);
    }

    std::vector<std::string> toLines() {
        return ParamOptions::toLines("SDInputOutput");
    }
};

class GenChainParamOptions : public ParamOptions {
public:
    GenChainParamOptions() {
        options_ = {
            {"InputCoordFile", ""},
            {"RMSDRefPDBFile", ""},
            {"SeqString", ""},
            {"RefPDBFile", ""},
            {"TypeInCoil", "GLY"},
            {"TypeInHelix", "LEU"},
            {"TypeInStrand", "VAL"}
        };
        requiredOptions_ = {
            "StartPDBFile",
            "SCTypeFrom",
            "SCCoordFrom"
        };
    }

    void setup(std::map<std::string, QA> qas) {
        std::string optname;
        std::string rawAnswer;

        // StartPDBFile
        ::setupPlainOption("StartPDBFile", qas, options_);

        // SCTypeFrom
        optname = "SCTypeFrom";
        rawAnswer = ::findRawAnswer(optname, qas);
        std::map<std::string, std::string> sctypemap {
            {"refpdb", "refpdb"},
            {"ss_type_in_refpdb", "sstypeinrefpdb"},
            {"userdefined", "seqstring"}
        };
        if (sctypemap.find(rawAnswer) == sctypemap.end()) {
            std::cerr << "Cease configuring! Invalid answer of "
                    << QAIndices[optname] << ": " << rawAnswer << std::endl;
            exit(1);
        }
        options_[optname] = sctypemap.at(rawAnswer);

        // SeqString
        ::setupPlainOption("SeqString", qas, options_);

        // SCCoordFrom
        ::setupPlainOption("SCCoordFrom", qas, options_);

        // RefPDBFile
        ::setupPlainOption("RefPDBFile", qas, options_);

        // InputCoordFile
        ::setupPlainOption("InputCoordFile", qas, options_);

        // RMSDRefPDBFile
        ::setupPlainOption("RMSDRefPDBFile", qas, options_);
    }

    std::vector<std::string> toLines() {
        return ParamOptions::toLines("GenChain");
    }
};

class SDParamOptions : public ParamOptions {
public:
    SDParamOptions() {
        options_ = {
            {"DoShake", "1"},
            {"TrajFile", ""},
            {"NTrajSteps", "0"},
            {"FixMainChain", "0"},
            {"NeighborListSteps", "50"},
            {"TimeStep", "0.002"},
            {"AnnealingScheme", ""}
        };
        requiredOptions_ = {
            "FrictionCoeff",
            "Temperatures",
        };
    }

    void setup(std::map<std::string, QA> qas) {
        // FrictionCoeff
        ::setupPlainOption("FrictionCoeff", qas, options_);

        std::string rawAnswer;
        std::string optname;

        // Temperatures, AnnealingScheme
        optname = "Temperatures";
        rawAnswer = ::findRawAnswer(optname, qas);
        if (!rawAnswer.empty()) {
            std::vector<std::string> fields;
            boost::split(fields, rawAnswer, boost::is_any_of(","));
            if (fields.size() == 1) {
                options_[optname] = fields[0];
            } else {
                options_[optname] = fields[1];
                std::string annscm = "0 " + fields[0] + " " + fields[1] + " " + fields[2] + " "
                                        + fields[3] + " " + fields[4];
                options_["AnnealingScheme"] = annscm;
            }
        }

        // TrajFile, NTrajSteps
        optname = "ZZZTrajOutput";
        rawAnswer = ::findRawAnswer(optname, qas);
        if (!rawAnswer.empty()) {
            std::vector<std::string> fields;
            boost::split(fields, rawAnswer, boost::is_any_of(","));
            if (fields.size() == 2) {
                options_["TrajFile"] = fields[0];
                options_["NTrajSteps"] = fields[1];
            }
        }
    }

    std::vector<std::string> toLines() {
        return ParamOptions::toLines("SD");
    }
};

class ForceFieldParamOptions : public ParamOptions {
public:
    ForceFieldParamOptions() {
        options_ = {
            {"CoreRegion", ""},
            {"CoreRegionBasePDB", ""},
            {"StrandPairRestrains", ""},
            {"EAttraction", "0"},
            {"KresStrandPair", "50"},
            {"LSWeight", "0.50"},
            {"PairPackingWeight", "0.32"},
            {"SideConfWeight", "2.40"},
            {"SCPackingWeight", "3.10"},
            {"LocalHBWeight", "0.60"},
            {"ShadowWeight", "0.00"},
            {"PhiPsiWeight", "2.00"},
            {"RAttractionOff", "100"},
            {"RAttractionSwitch", "80"},
            {"StericWeight", "1.40"},
            {"SideChainMMSteric", "0"},
            {"SubStructureRestrains", ""},
            {"RgRestrain", "5.0 -300.0"},
            {"TotalRMSDRestrain", ""},
            {"RMSDRestrainedRegion", ""},
            {"DisRestrains", ""},
            {"ExtendRestrains", ""},
            {"HelixRestrains", ""},
            {"LoopRegionWeights", ""}
        };
        requiredOptions_ = {
        };
    }

    void setup(std::map<std::string, QA> qas) {
        std::string optname;
        std::string rawAnswer;
        std::string optstr = "";

        // LocalHBWeight
        optname = "LocalHBWeight";
        rawAnswer = ::findRawAnswer(optname, qas);
        if (rawAnswer == "no") {
            options_[optname] = "0.00";
        }

        // SideConfWeight, SCPackingWeight
        optname = "ZZZSCEnergyTermsWeight";
        rawAnswer = ::findRawAnswer(optname, qas);
        if (rawAnswer == "no") {
            options_["SideConfWeight"] = "0.24";
            options_["SCPackingWeight"] = "0.00";
        }

        // RgRestrain
        optname = "RgRestrain";
        rawAnswer = ::findRawAnswer(optname, qas);
        if (rawAnswer == "no") {
            options_[optname] = "1000.0 0.0";
        } else {
            options_[optname] = "5.0 -300.0";
            // ZZZRgRestrainRegionFrom
            optname = "ZZZRgRestrainRegionFrom";
            rawAnswer = ::findRawAnswer(optname, qas);
            if ("ss_in_refpdb" == rawAnswer) {
                // CoreRegionBasePDB
                optname = "CoreRegionBasePDB";
                std::string corepdbfile = ::findRawAnswer("RefPDBFile", qas);
                if (corepdbfile.empty()) {
                    corepdbfile = ::findRawAnswer("StartPDBFile", qas);
                }
                options_[optname] = corepdbfile;
            } else if ("userdefined" == rawAnswer) {
                // CoreRegion
                optname = "CoreRegion";
                rawAnswer = ::findRawAnswer(optname, qas);
                std::vector<std::string> fields;
                boost::split(fields, rawAnswer, boost::is_any_of(","));
                std::string optstr = "";
                for (auto & field : fields) {
                    std::vector<std::string> numbers;
                    boost::split(numbers, field, boost::is_any_of("-"));
                    int start = 0;
                    int end = 0;
                    if (1 == numbers.size()) {
                        int site = std::stoi(numbers[0]);
                        start = site - 1; // sequence number to index
                        end = start; // single position, same as start
                    } else {
                        int site1 = std::stoi(numbers[0]);
                        start = site1 - 1; // sequence number to index
                        int site2 = std::stoi(numbers[1]);
                        end = site2 - 1; // sequence number to index
                    }
                    optstr.append("0 "); // only single chain at the moment
                    optstr.append(std::to_string(start) + " ");
                    optstr.append(std::to_string(end) + " ");
                }
                boost::algorithm::trim(optstr);
                options_[optname] = optstr;
            } else {
                options_["CoreRegion"] = "";
                options_["CoreRegionBasePDB"] = "";
            }
        }

        // TotalRMSDRestrain
        optname = "TotalRMSDRestrain";
        rawAnswer = ::findRawAnswer(optname, qas);
        if ("none" == rawAnswer) {
            options_[optname] = "";
            options_["RMSDRestrainedRegion"] = "";
        } else {
            std::string rmsd0 = ::findRawAnswer("ZZZRMSDRestrainValue", qas);
            std::string krmsdrst = "1000.0";
            std::map<std::string, std::string> rmsdrsttypemap {
                {"backbone_only", "backboneonly"},
                {"sse_backbone_only", "sseonly"},
                {"all_atoms", "allatoms"}
            };
            if (rmsdrsttypemap.find(rawAnswer) == rmsdrsttypemap.end()) {
                std::cerr << "Cease configuring! Invalid answer of "
                        << QAIndices[optname] << ": " << rawAnswer << std::endl;
                exit(1);
            }
            std::string rmsdrsttype = rmsdrsttypemap.at(rawAnswer);
            optstr = "total_mode " + rmsd0 + " " + krmsdrst + " " + rmsdrsttype;
            options_[optname] = optstr;

            // RMSDRestrainedRegion
            optname = "RMSDRestrainedRegion";
            optstr = "";
            rawAnswer = ::findRawAnswer(optname, qas);
            if (rawAnswer.empty()) {
                options_[optname] = "";
            } else {
                std::vector<std::string> fields;
                boost::split(fields, rawAnswer, boost::is_any_of(","));
                for (auto & field : fields) {
                    std::vector<std::string> numbers;
                    boost::split(numbers, field, boost::is_any_of("-"));
                    int start = 0;
                    int end = 0;
                    if (1 == numbers.size()) {
                        int site = std::stoi(numbers[0]);
                        start = site - 1; // sequence number to index
                        end = start; // single position, same as start
                    } else {
                        int site1 = std::stoi(numbers[0]);
                        int site2 = std::stoi(numbers[1]);
                        start = site1 - 1; // sequence number to index
                        end = site2 - 1; // sequence number to index
                    }
                    optstr = optstr.append("0 "); // only single chain at the moment
                    optstr = optstr.append(std::to_string(start) + " ");
                    optstr = optstr.append(std::to_string(end) + " ");
                }
                boost::algorithm::trim(optstr);
            }
            options_[optname] = optstr;
        }

        // HelixRestrains
        optname = "HelixRestrains";
        optstr = "";
        rawAnswer = ::findRawAnswer(optname, qas);
        if (rawAnswer.empty()) {
            options_[optname] = "";
        } else {
            std::string rstparam = "0.99 1000.0";
            std::vector<std::string> fields;
            boost::split(fields, rawAnswer, boost::is_any_of(","));
            for (auto & field : fields) {
                std::vector<std::string> numbers;
                boost::split(numbers, field, boost::is_any_of("-"));
                int start = 0;
                int end = 0;
                if (1 == numbers.size()) {
                    int site = std::stoi(numbers[0]);
                    start = site - 1; // sequence number to index
                    end = start; // single position, same as start
                } else {
                    int site1 = std::stoi(numbers[0]);
                    int site2 = std::stoi(numbers[1]);
                    start = site1 - 1; // sequence number to index
                    end = site2 - 1; // sequence number to index
                }
                optstr = optstr.append("0 "); // only single chain at the moment
                optstr = optstr.append(std::to_string(start) + " ");
                optstr = optstr.append(std::to_string(end) + " ");
                optstr = optstr.append(rstparam + " ");
            }
            boost::algorithm::trim(optstr);
        }
        options_[optname] = optstr;

        // ExtendRestrains
        optname = "ExtendRestrains";
        optstr = "";
        rawAnswer = ::findRawAnswer(optname, qas);
        if (rawAnswer.empty()) {
            options_[optname] = "";
        } else {
            std::vector<std::string> fields;
            boost::split(fields, rawAnswer, boost::is_any_of(","));
            for (auto & field : fields) {
                std::vector<std::string> numbers;
                boost::split(numbers, field, boost::is_any_of("-"));
                int start = 0;
                int end = 0;
                if (1 == numbers.size()) {
                    int site = std::stoi(numbers[0]);
                    start = site - 1; // sequence number to index
                    end = start; // single position, same as start
                } else {
                    int site1 = std::stoi(numbers[0]);
                    int site2 = std::stoi(numbers[1]);
                    start = site1 - 1; // sequence number to index
                    end = site2 - 1; // sequence number to index
                }
                optstr = optstr.append("0 "); // only single chain at the moment
                optstr = optstr.append(std::to_string(start) + " ");
                optstr = optstr.append(std::to_string(end) + " ");
            }
            boost::algorithm::trim(optstr);
        }
        options_[optname] = optstr;

        // LoopRegionWeights
        optname = "LoopRegionWeights";
        optstr = "";
        rawAnswer = ::findRawAnswer(optname, qas);
        if (rawAnswer.empty()) {
            options_[optname] = "";
        } else {
            std::vector<std::string> fields;
            boost::split(fields, rawAnswer, boost::is_any_of(","));
            std::string rstweight = "0.001";
            for (auto & field : fields) {
                std::vector<std::string> numbers;
                boost::split(numbers, field, boost::is_any_of("-"));
                int start = 0;
                int end = 0;
                if (1 == numbers.size()) {
                    int site = std::stoi(numbers[0]);
                    start = site - 1; // sequence number to index
                    end = start; // single position, same as start
                } else {
                    int site1 = std::stoi(numbers[0]);
                    int site2 = std::stoi(numbers[1]);
                    start = site1 - 1; // sequence number to index
                    end = site2 - 1; // sequence number to index
                }
                optstr = optstr.append("0 "); // only single chain at the moment
                optstr = optstr.append(std::to_string(start) + " ");
                optstr = optstr.append(std::to_string(end) + " ");
                optstr = optstr.append(rstweight + " ");
            }
            boost::algorithm::trim(optstr);
        }
        options_[optname] = optstr;
    }

    std::vector<std::string> toLines() {
        return ParamOptions::toLines("ForceField");
    }
};


void configure(const std::string & inputConfigFile,
              const std::string & outputParamFile)
{
    std::map<std::string, QA> qas = parseQALines(inputConfigFile);

    SDIOParamOptions sdiopo;
    SDParamOptions sdpo;
    GenChainParamOptions gcpo;
    ForceFieldParamOptions ffpo;

    sdiopo.setup(qas);
    sdpo.setup(qas);
    gcpo.setup(qas);
    ffpo.setup(qas);

    bool setupGood = sdiopo.isSufficient() && sdpo.isSufficient()
                        && gcpo.isSufficient() && ffpo.isSufficient();
    if (! setupGood) {
        std::cerr << "Problem(s) detected in the configuration file, please check it." << std::endl;
        exit(1);
    }

    std::vector<std::string> paramLines;
    std::vector<std::string> lines;

    lines = sdiopo.toLines();
    for (auto & line : lines) {
        paramLines.push_back(line);
    }
    lines = sdpo.toLines();
    for (auto & line : lines) {
        paramLines.push_back(line);
    }
    lines = gcpo.toLines();
    for (auto & line : lines) {
        paramLines.push_back(line);
    }
    lines = ffpo.toLines();
    for (auto & line : lines) {
        paramLines.push_back(line);
    }

    std::ofstream ofs(outputParamFile);
    if (!ofs.is_open()) {
        std::cerr << "Failed to open the output file <"
                << outputParamFile << ">." << std::endl;
        exit(1);
    }
    for (auto & line : paramLines) {
        ofs << line << std::endl;
    }
    ofs.close();
    std::cout << "The parameters have been written to file <"
            << outputParamFile << "> successfully." << std::endl;
}

void printUsage(const std::string & selfname)
{
    std::cout << "Usage:" << std::endl;
    std::cout << "$ " << selfname << " <inputConfigFile> <outputParamFile>"
            << std::endl;
}

int main(int argc, char** argv)
{
    std::string selfname = argv[0];

    if (argc != 3) {
        std::cerr << "Wrong number of arguments." << std::endl;
        printUsage(selfname);
        return 1;
    }

    std::string inputConfigFile(argv[1]);
    std::string outputParamFile(argv[2]);

    configure(inputConfigFile, outputParamFile);

    return 0;
}
