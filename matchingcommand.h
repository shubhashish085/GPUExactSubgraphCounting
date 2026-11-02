#ifndef GSC_MATCHINGCOMMAND_H
#define GSC_MATCHINGCOMMAND_H

#include "commandparser.h"
#include <map>
#include <iostream>
enum OptionKeyword {
    QueryGraphFile = 0,     // -q, The query graph file path, compulsive parameter
    DataGraphFile = 1,      // -d, The data graph file path, compulsive parameter
    OutputFilePath = 2      // -output, the output file path, compulsive parameter
};

class MatchingCommand : public CommandParser{
private:
    std::map<OptionKeyword, std::string> options_key;
    std::map<OptionKeyword, std::string> options_value;

private:
    void processOptions();

public:
    MatchingCommand(int argc, char **argv);

    std::string getDataGraphFilePath() {
        return options_value[OptionKeyword::DataGraphFile];
    }

    std::string getQueryGraphFilePath() {
        return options_value[OptionKeyword::QueryGraphFile];
    }

    std::string getOutputFilePath(){
        return options_value[OptionKeyword::OutputFilePath];
    }
};



#endif