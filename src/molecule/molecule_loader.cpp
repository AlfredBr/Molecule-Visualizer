#include "molecule/molecule_loader.h"
#include "molecule/molecule_db.h"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#include <limits.h>
#endif

namespace fs = std::filesystem;

namespace {

struct Json {
    enum Type { Null, Bool, Number, String, Array, Object } type = Null;
    bool boolean = false;
    double number = 0.0;
    std::string string;
    std::vector<Json> array;
    std::map<std::string, Json> object;

    const Json* get(const char* key) const {
        if (type != Object) return nullptr;
        auto it = object.find(key);
        return it == object.end() ? nullptr : &it->second;
    }
};

class JsonParser {
public:
    explicit JsonParser(const std::string& text) : text_(text) {}

    Json parse() {
        Json result = value();
        whitespace();
        if (pos_ != text_.size()) fail("unexpected trailing content");
        return result;
    }

private:
    const std::string& text_;
    size_t pos_ = 0;

    [[noreturn]] void fail(const std::string& message) const {
        throw std::runtime_error(message + " at byte " + std::to_string(pos_));
    }

    void whitespace() {
        while (pos_ < text_.size() && (text_[pos_] == ' ' || text_[pos_] == '\t' ||
               text_[pos_] == '\r' || text_[pos_] == '\n')) ++pos_;
    }

    bool consume(char c) {
        whitespace();
        if (pos_ < text_.size() && text_[pos_] == c) { ++pos_; return true; }
        return false;
    }

    Json value() {
        whitespace();
        if (pos_ >= text_.size()) fail("expected JSON value");
        char c = text_[pos_];
        if (c == '{') return object();
        if (c == '[') return array();
        if (c == '"') { Json j; j.type = Json::String; j.string = string(); return j; }
        if (c == 't') return literal("true", Json::Bool, true);
        if (c == 'f') return literal("false", Json::Bool, false);
        if (c == 'n') return literal("null", Json::Null, false);
        if (c == '-' || (c >= '0' && c <= '9')) return number();
        fail("unexpected character");
    }

    Json literal(const char* word, Json::Type type, bool boolean) {
        size_t length = std::char_traits<char>::length(word);
        if (text_.compare(pos_, length, word) != 0) fail("invalid literal");
        pos_ += length;
        Json j; j.type = type; j.boolean = boolean; return j;
    }

    Json number() {
        const char* begin = text_.c_str() + pos_;
        char* end = nullptr;
        errno = 0;
        double result = std::strtod(begin, &end);
        if (end == begin || errno == ERANGE || !std::isfinite(result)) fail("invalid number");
        pos_ += static_cast<size_t>(end - begin);
        Json j; j.type = Json::Number; j.number = result; return j;
    }

    std::string string() {
        if (!consume('"')) fail("expected string");
        std::string out;
        while (pos_ < text_.size()) {
            unsigned char c = static_cast<unsigned char>(text_[pos_++]);
            if (c == '"') return out;
            if (c < 0x20) fail("control character in string");
            if (c != '\\') { out.push_back(static_cast<char>(c)); continue; }
            if (pos_ >= text_.size()) fail("unterminated escape");
            char e = text_[pos_++];
            switch (e) {
                case '"': out.push_back('"'); break; case '\\': out.push_back('\\'); break;
                case '/': out.push_back('/'); break; case 'b': out.push_back('\b'); break;
                case 'f': out.push_back('\f'); break; case 'n': out.push_back('\n'); break;
                case 'r': out.push_back('\r'); break; case 't': out.push_back('\t'); break;
                case 'u': {
                    if (pos_ + 4 > text_.size()) fail("short unicode escape");
                    unsigned code = 0;
                    for (int i = 0; i < 4; ++i) {
                        char h = text_[pos_++];
                        code <<= 4;
                        if (h >= '0' && h <= '9') code += h - '0';
                        else if (h >= 'a' && h <= 'f') code += h - 'a' + 10;
                        else if (h >= 'A' && h <= 'F') code += h - 'A' + 10;
                        else fail("invalid unicode escape");
                    }
                    if (code <= 0x7f) out.push_back(static_cast<char>(code));
                    else if (code <= 0x7ff) {
                        out.push_back(static_cast<char>(0xc0 | (code >> 6)));
                        out.push_back(static_cast<char>(0x80 | (code & 0x3f)));
                    } else {
                        out.push_back(static_cast<char>(0xe0 | (code >> 12)));
                        out.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3f)));
                        out.push_back(static_cast<char>(0x80 | (code & 0x3f)));
                    }
                    break;
                }
                default: fail("invalid escape");
            }
        }
        fail("unterminated string");
    }

    Json array() {
        consume('['); Json j; j.type = Json::Array;
        if (consume(']')) return j;
        do { j.array.push_back(value()); } while (consume(','));
        if (!consume(']')) fail("expected ']'");
        return j;
    }

    Json object() {
        consume('{'); Json j; j.type = Json::Object;
        if (consume('}')) return j;
        do {
            whitespace();
            if (pos_ >= text_.size() || text_[pos_] != '"') fail("expected object key");
            std::string key = string();
            if (!consume(':')) fail("expected ':'");
            if (!j.object.emplace(key, value()).second) fail("duplicate object key");
        } while (consume(','));
        if (!consume('}')) fail("expected '}'");
        return j;
    }
};

bool requiredString(const Json& root, const char* key, std::string* out, std::string* error) {
    const Json* value = root.get(key);
    if (!value || value->type != Json::String || value->string.empty()) {
        *error = std::string("'") + key + "' must be a nonempty string"; return false;
    }
    *out = value->string; return true;
}

bool optionalString(const Json& root, const char* key, std::string* out, std::string* error) {
    const Json* value = root.get(key);
    if (!value) return true;
    if (value->type != Json::String) { *error = std::string("'") + key + "' must be a string"; return false; }
    *out = value->string; return true;
}

int atomType(const std::string& symbol) {
    static const char* symbols[ATOM_TYPE_COUNT] = {
        "H","C","N","O","P","S","Cl","Br","F","I","Na","Si","B","Fe","Cu",
        "Al","Ti","Pt","Re","Xe","Mo","W","Se","Ge","Ga"
    };
    for (int i = 0; i < ATOM_TYPE_COUNT; ++i) if (symbol == symbols[i]) return i;
    return -1;
}

float atomRadius(int type) {
    static const float radii[ATOM_TYPE_COUNT] = {
        .25f,.40f,.38f,.35f,.45f,.45f,.45f,.50f,.35f,.55f,.55f,.48f,.42f,
        .55f,.50f,.50f,.52f,.55f,.62f,.62f,.55f,.58f,.50f,.50f,.50f
    };
    return radii[type];
}

int categoryFromName(const std::string& name) {
    for (int i = 0; i < CAT_COUNT; ++i) if (name == molecule_get_category_name(i)) return i;
    return -1;
}

bool validId(const std::string& id) {
    if (id.empty()) return false;
    for (char c : id) if (!((c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '-')) return false;
    return true;
}

void center(Molecule* molecule) {
    if (molecule->numAtoms == 0) return;
    float x = 0, y = 0, z = 0;
    for (int i = 0; i < molecule->numAtoms; ++i) {
        x += molecule->atoms[i].x; y += molecule->atoms[i].y; z += molecule->atoms[i].z;
    }
    x /= molecule->numAtoms; y /= molecule->numAtoms; z /= molecule->numAtoms;
    for (int i = 0; i < molecule->numAtoms; ++i) {
        molecule->atoms[i].x -= x; molecule->atoms[i].y -= y; molecule->atoms[i].z -= z;
    }
}

} // namespace

bool molecule_parse_json_file(const std::string& path, ExternalMoleculeRecord* record,
                              std::string* error) {
    try {
        std::ifstream input(path, std::ios::binary);
        if (!input) { *error = "could not open file"; return false; }
        std::ostringstream buffer; buffer << input.rdbuf();
        Json root = JsonParser(buffer.str()).parse();
        if (root.type != Json::Object) { *error = "root must be an object"; return false; }
        const Json* version = root.get("schemaVersion");
        if (!version || version->type != Json::Number || version->number != 1) {
            *error = "'schemaVersion' must be 1"; return false;
        }
        ExternalMoleculeRecord result;
        if (!requiredString(root, "id", &result.id, error) ||
            !requiredString(root, "name", &result.name, error) ||
            !optionalString(root, "formula", &result.formula, error) ||
            !optionalString(root, "description", &result.description, error) ||
            !optionalString(root, "longDescription", &result.longDescription, error) ||
            !optionalString(root, "structureType", &result.structureType, error)) return false;
        if (!validId(result.id)) { *error = "'id' may contain only lowercase letters, digits, and hyphens"; return false; }
        if (result.name.size() >= sizeof(result.molecule.name)) { *error = "'name' is too long"; return false; }
        if (result.formula.size() >= sizeof(result.molecule.formula)) { *error = "'formula' is too long"; return false; }
        std::string category;
        if (!requiredString(root, "category", &category, error)) return false;
        result.category = categoryFromName(category);
        if (result.category < 0) { *error = "unknown category '" + category + "'"; return false; }
        if (const Json* over = root.get("override")) {
            if (over->type != Json::Bool) { *error = "'override' must be a boolean"; return false; }
            result.overrideExisting = over->boolean;
        }
        const Json* atoms = root.get("atoms");
        if (!atoms || atoms->type != Json::Array || atoms->array.empty()) { *error = "'atoms' must be a nonempty array"; return false; }
        if (atoms->array.size() > MAX_ATOMS) { *error = "atom count exceeds MAX_ATOMS"; return false; }
        for (const Json& atom : atoms->array) {
            if (atom.type != Json::Object) { *error = "each atom must be an object"; return false; }
            std::string symbol;
            if (!requiredString(atom, "element", &symbol, error)) return false;
            int type = atomType(symbol);
            if (type < 0) { *error = "unsupported element '" + symbol + "'"; return false; }
            const Json* position = atom.get("position");
            if (!position || position->type != Json::Array || position->array.size() != 3) { *error = "atom position must contain three numbers"; return false; }
            Atom& output = result.molecule.atoms[result.molecule.numAtoms++];
            float* coordinates[3] = { &output.x, &output.y, &output.z };
            for (int i = 0; i < 3; ++i) {
                if (position->array[i].type != Json::Number || std::fabs(position->array[i].number) > 10000.0) { *error = "atom coordinates must be finite and within 10000 angstroms"; return false; }
                *coordinates[i] = static_cast<float>(position->array[i].number);
            }
            output.type = type; output.radius = atomRadius(type);
        }
        const Json* bonds = root.get("bonds");
        std::set<std::pair<int,int>> seen;
        if (bonds) {
            if (bonds->type != Json::Array || bonds->array.size() > MAX_BONDS) { *error = "'bonds' must be an array within MAX_BONDS"; return false; }
            for (const Json& bond : bonds->array) {
                if (bond.type != Json::Object) { *error = "each bond must be an object"; return false; }
                const Json* pair = bond.get("atoms"); const Json* order = bond.get("order");
                if (!pair || pair->type != Json::Array || pair->array.size() != 2 || !order || order->type != Json::Number) { *error = "bond requires two atom indices and an order"; return false; }
                if (pair->array[0].type != Json::Number || pair->array[1].type != Json::Number) { *error = "bond indices must be integers"; return false; }
                int a = static_cast<int>(pair->array[0].number), b = static_cast<int>(pair->array[1].number), o = static_cast<int>(order->number);
                if (a != pair->array[0].number || b != pair->array[1].number || o != order->number) { *error = "bond indices and order must be integers"; return false; }
                if (a < 0 || b < 0 || a >= result.molecule.numAtoms || b >= result.molecule.numAtoms || a == b) { *error = "bond contains an invalid atom index"; return false; }
                if (o < 1 || o > 3) { *error = "bond order must be 1, 2, or 3"; return false; }
                auto key = std::minmax(a, b);
                if (!seen.emplace(key.first, key.second).second) { *error = "duplicate bond"; return false; }
                result.molecule.bonds[result.molecule.numBonds++] = {a, b, o};
            }
        }
        std::snprintf(result.molecule.name, sizeof(result.molecule.name), "%s", result.name.c_str());
        std::snprintf(result.molecule.formula, sizeof(result.molecule.formula), "%s", result.formula.c_str());
        center(&result.molecule);
        result.source = path;
        *record = std::move(result);
        return true;
    } catch (const std::exception& ex) {
        *error = ex.what(); return false;
    }
}

MoleculeLoadResult molecule_load_json_directories(const std::vector<std::string>& directories) {
    MoleculeLoadResult result;
    for (const std::string& directory : directories) {
        std::error_code ec;
        if (!fs::is_directory(fs::u8path(directory), ec)) continue;
        std::vector<fs::path> files;
        for (const auto& entry : fs::directory_iterator(fs::u8path(directory), ec))
            if (entry.is_regular_file() && entry.path().extension() == ".json") files.push_back(entry.path());
        std::sort(files.begin(), files.end());
        for (const fs::path& file : files) {
            ExternalMoleculeRecord record; std::string error;
            std::string path = file.u8string();
            if (molecule_parse_json_file(path, &record, &error)) result.records.push_back(std::move(record));
            else result.errors.push_back(path + ": " + error);
        }
    }
    return result;
}

std::string molecule_application_directory() {
#ifdef _WIN32
    wchar_t buffer[32768]; DWORD length = GetModuleFileNameW(nullptr, buffer, 32768);
    if (length > 0 && length < 32768) return fs::path(buffer).parent_path().u8string();
#elif defined(__APPLE__)
    char buffer[PATH_MAX]; uint32_t size = sizeof(buffer);
    if (_NSGetExecutablePath(buffer, &size) == 0) return fs::weakly_canonical(buffer).parent_path().u8string();
#endif
    return fs::current_path().u8string();
}

std::string molecule_user_directory() {
#ifdef _WIN32
    const char* local = std::getenv("LOCALAPPDATA");
    return (local ? fs::u8path(local) : fs::current_path()) .append("MolVis").append("molecules").u8string();
#else
    const char* home = std::getenv("HOME");
    fs::path base = home ? fs::u8path(home) : fs::current_path();
#ifdef __APPLE__
    return base.append("Library").append("Application Support").append("MolVis").append("molecules").u8string();
#else
    return base.append(".config").append("molvis").append("molecules").u8string();
#endif
#endif
}

bool molecule_ensure_user_directory(std::string* error) {
    std::error_code ec; fs::create_directories(fs::u8path(molecule_user_directory()), ec);
    if (ec) { if (error) *error = ec.message(); return false; }
    return true;
}
