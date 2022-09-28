//
// Programmer: Created by Alex Stauffer on 07.09.21.
// File: SidebandUtils.cpp
// Purpose: Tools to create help SidebandFit scripts import the data
//

#include "SidebandUtils.h"
// #include "LatexTableWriter.h"

// This function takes as input a terminal command and return the standard output
std::vector<std::string> GetStdoutFromCommand(std::string cmd) {

  std::vector<std::string> data;
  FILE * stream;
  const int max_buffer = 256;
  char buffer[max_buffer];
  cmd.append(" 2>&1");

  stream = popen(cmd.c_str(), "r");
  if (stream) {
    while (!feof(stream))
      if (fgets(buffer, max_buffer, stream) != NULL) data.push_back(buffer);
    pclose(stream);
  }
  return data;
}

// Extracts and combines data from a RooWorkspace from the given Satoshi Nrs and input path from a given S1 bin
std::unique_ptr<RooDataSet> GetDataFromWorkspaces(std::string Path, std::string DataType, int MinRunNr, int MaxRunNr, int nr_S1bin) {

  // Extract list of filenames with/out source (datatype)
  std::string command = "ls " + Path + "WarioRun*" + DataType + ".root";
  std::vector <std::string> pre_file_list = GetStdoutFromCommand(command);
  int run_ID_string_length = 4;
  int start = InputPath.length() + 8;  // getting the starting index of run number 00xx

  // Looping through the files and keeping the ones with appropriate run number
  // mainly to get the number of files to be processed for progress information...
  std::vector <std::string> file_list;
  for (std::vector<std::string>::iterator t = pre_file_list.begin(); t != pre_file_list.end(); ++t) {
    const int run_nr_temp = std::stoi(t->substr(start, run_ID_string_length));
    if (run_nr_temp >= MinRunNr && run_nr_temp <= MaxRunNr) {
      file_list.push_back(t->substr(0, t->length() - 1));
    }
  }
  int files_to_be_processed = file_list.size();

  // Pointer to the combined data, we use a smart pointer here so we
  // don't have to remember to call `delete`.
  std::unique_ptr<RooDataSet> Data_in_S1Bin;

  // Looping through all files
  int processed_files = 0;
  for (std::vector<std::string>::iterator name = file_list.begin(); name != file_list.end(); ++name) {
    // The TFile can be created on the stack
    TFile file(name->c_str());
    processed_files += 1;

    // We have to soon delete the RooWorkspace that is retrieved from the file,
    // so let's wrap it into a unique_ptr such that this happens automatically.
    std::unique_ptr<RooWorkspace> workspace{file.Get<RooWorkspace>("w")};

    // Get a non-owning pointer to the dataset (don't delete this,
    // it is owned by the workspace).
    auto * DataToAdd = static_cast<RooDataSet*>(workspace->data(("S1bin_" + to_string(nr_S1bin) + "_single").c_str()));

    if(processed_files == 1) {
      // If this is the first file, we create our combined dataset by
      // cloning the dataset in the first file.
      // std::cout << "Initialising:" << name->c_str() << " (1/" << files_to_be_processed << ")" << std::endl;
      Data_in_S1Bin.reset(static_cast<RooDataSet*>(workspace->data(("S1bin_" + to_string(nr_S1bin) + "_single").c_str())->Clone()));
    } else {
      // Otherwise, just append the data.
      // std::cout << "Adding:" << name->c_str() << " (" << processed_files << "/" << files_to_be_processed << ")" << std::endl;
      Data_in_S1Bin->append(*DataToAdd);
    }
  }
  // Changing the name of the DataSet to incorporate the datatype
  Data_in_S1Bin->SetName(("S1bin_" + to_string(nr_S1bin) + "_" + DataType).c_str());
  // Returning the unique_ptr to the Dataset
  return Data_in_S1Bin;
}
