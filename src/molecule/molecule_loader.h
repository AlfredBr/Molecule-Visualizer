#ifndef MOLECULE_LOADER_H
#define MOLECULE_LOADER_H

#include "molecule/molecule_types.h"
#include <string>
#include <vector>

struct ExternalMoleculeRecord {
    std::string id;
    std::string name;
    std::string formula;
    int category = 0;
    std::string description;
    std::string longDescription;
    std::string structureType;
    std::string source;
    bool overrideExisting = false;
    Molecule molecule = {};
};

struct MoleculeLoadResult {
    std::vector<ExternalMoleculeRecord> records;
    std::vector<std::string> errors;
};

MoleculeLoadResult molecule_load_json_directories(const std::vector<std::string>& directories);
bool molecule_parse_json_file(const std::string& path, ExternalMoleculeRecord* record,
                              std::string* error);
std::string molecule_application_directory();
std::string molecule_user_directory();
bool molecule_ensure_user_directory(std::string* error);

#endif
