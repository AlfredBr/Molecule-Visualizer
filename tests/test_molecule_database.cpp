#include "molecule/molecule_db.h"
#include "molecule/molecule_loader.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <cstring>

namespace fs = std::filesystem;

static int failures = 0;
#define CHECK(condition) do { if (!(condition)) { std::cerr << "FAIL " << __FILE__ << ':' << __LINE__ << ": " #condition "\n"; ++failures; } } while (0)

int main() {
    CHECK(molecule_get_count() > 290);
    std::set<std::string> stableIds;
    for (int index = 0; index < molecule_get_count(); ++index) {
        CHECK(stableIds.insert(molecule_get_id(index)).second);
        Molecule molecule = {};
        molecule_build(index, &molecule);
        CHECK(std::memchr(molecule.name, '\0', sizeof(molecule.name)) != nullptr);
        CHECK(std::memchr(molecule.formula, '\0', sizeof(molecule.formula)) != nullptr);
        CHECK(molecule.numAtoms >= 0 && molecule.numAtoms <= MAX_ATOMS);
        CHECK(molecule.numBonds >= 0 && molecule.numBonds <= MAX_BONDS);
        for (int atom = 0; atom < molecule.numAtoms; ++atom)
            CHECK(molecule.atoms[atom].type >= 0 && molecule.atoms[atom].type < ATOM_TYPE_COUNT);
        for (int bond = 0; bond < molecule.numBonds; ++bond) {
            CHECK(molecule.bonds[bond].atom1 >= 0 && molecule.bonds[bond].atom1 < molecule.numAtoms);
            CHECK(molecule.bonds[bond].atom2 >= 0 && molecule.bonds[bond].atom2 < molecule.numAtoms);
            CHECK(molecule.bonds[bond].order >= 1 && molecule.bonds[bond].order <= 3);
        }
    }
    const char* ids[] = { "molybdenum-disulfide", "tungsten-disulfide", "tungsten-diselenide",
                          "crystalline-silicon", "silicon-germanium", "silicon-carbide-3c",
                          "gallium-nitride-wurtzite" };
    for (const char* id : ids) {
        int index = molecule_find_by_id(id);
        CHECK(index >= 0);
        if (index >= 0) {
            Molecule molecule = {};
            molecule_build(index, &molecule);
            CHECK(molecule.numAtoms > 0);
            CHECK(molecule.numAtoms <= MAX_ATOMS);
            CHECK(molecule.numBonds <= MAX_BONDS);
            for (int i = 0; i < molecule.numBonds; ++i) {
                CHECK(molecule.bonds[i].atom1 >= 0 && molecule.bonds[i].atom1 < molecule.numAtoms);
                CHECK(molecule.bonds[i].atom2 >= 0 && molecule.bonds[i].atom2 < molecule.numAtoms);
            }
        }
    }

    fs::path directory = fs::temp_directory_path() / "molvis-database-test";
    fs::create_directories(directory);
    fs::path valid = directory / "valid.json";
    {
        std::ofstream output(valid);
        output << R"({
          "schemaVersion": 1,
          "id": "test-water",
          "name": "Test Water",
          "formula": "H2O",
          "category": "Simple Molecules",
          "atoms": [
            {"element":"O", "position":[0,0,0]},
            {"element":"H", "position":[0.76,0.59,0]},
            {"element":"H", "position":[-0.76,0.59,0]}
          ],
          "bonds": [
            {"atoms":[0,1], "order":1},
            {"atoms":[0,2], "order":1}
          ]
        })";
    }
    ExternalMoleculeRecord record; std::string error;
    CHECK(molecule_parse_json_file(valid.u8string(), &record, &error));
    CHECK(record.id == "test-water");
    CHECK(record.molecule.numAtoms == 3);
    CHECK(record.molecule.numBonds == 2);

    fs::path invalid = directory / "invalid.json";
    { std::ofstream output(invalid); output << R"({"schemaVersion":1,"id":"Bad ID"})"; }
    CHECK(!molecule_parse_json_file(invalid.u8string(), &record, &error));
    CHECK(!error.empty());
    fs::remove_all(directory);

    if (failures == 0) std::cout << "All molecule database tests passed\n";
    return failures == 0 ? 0 : 1;
}
