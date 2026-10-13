#include "chem/kumo/fchk.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <zlib.h>

namespace chem::kumo {
namespace {

std::string trim_right(std::string s) {
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r')) {
        s.pop_back();
    }
    return s;
}

std::string normalize_fortran_exponent(std::string token) {
    for (auto& c : token) {
        if (c == 'D' || c == 'd') c = 'E';
    }
    return token;
}

long long parse_integer_token(const std::string& token, const std::string& label) {
    std::size_t consumed = 0;
    long long value = 0;
    try {
        value = std::stoll(token, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid integer token in fchk record " + label + ": " + token);
    }
    if (consumed != token.size()) {
        throw std::runtime_error("trailing characters in integer token in fchk record " + label + ": " + token);
    }
    return value;
}

double parse_real_token(const std::string& raw_token, const std::string& label) {
    const auto token = normalize_fortran_exponent(raw_token);
    std::size_t consumed = 0;
    double value = 0.0;
    try {
        value = std::stod(token, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid real token in fchk record " + label + ": " + raw_token);
    }
    if (consumed != token.size() || !std::isfinite(value)) {
        throw std::runtime_error("invalid real token in fchk record " + label + ": " + raw_token);
    }
    return value;
}

template <typename Parse>
void read_array_tokens(std::istream& input,
                       std::size_t count,
                       const std::string& label,
                       Parse&& parse) {
    for (std::size_t i = 0; i < count; ++i) {
        std::string token;
        if (!(input >> token)) {
            throw std::runtime_error("truncated array in fchk: " + label);
        }
        parse(token);
    }
    std::string remainder;
    std::getline(input, remainder);
    if (remainder.find_first_not_of(" \t\r") != std::string::npos) {
        throw std::runtime_error("too many values on final array line in fchk: " + label);
    }
}

// Inflate a complete gzip stream. Z_STREAM_END is mandatory: reaching the
// physical end of the file is not sufficient to validate CRC32 and ISIZE.
std::string decompress_gzip(const std::filesystem::path& path) {
    std::ifstream source(path, std::ios::binary);
    if (!source) throw std::runtime_error("cannot open gzip formatted checkpoint: " + path.string());

    z_stream stream{};
    if (inflateInit2(&stream, MAX_WBITS + 16) != Z_OK) {
        throw std::runtime_error("cannot initialize gzip decoder: " + path.string());
    }
    std::array<unsigned char, 65536> input_buffer{};
    std::array<unsigned char, 65536> output_buffer{};
    std::string contents;
    try {
        bool complete = false;
        while (!complete) {
            if (stream.avail_in == 0) {
                source.read(reinterpret_cast<char*>(input_buffer.data()),
                            static_cast<std::streamsize>(input_buffer.size()));
                const auto count = source.gcount();
                if (count == 0) {
                    throw std::runtime_error("truncated gzip formatted checkpoint: " + path.string());
                }
                stream.next_in = input_buffer.data();
                stream.avail_in = static_cast<uInt>(count);
            }
            stream.next_out = output_buffer.data();
            stream.avail_out = static_cast<uInt>(output_buffer.size());
            const int status = inflate(&stream, Z_NO_FLUSH);
            const auto produced = output_buffer.size() - stream.avail_out;
            contents.append(reinterpret_cast<const char*>(output_buffer.data()), produced);
            if (status == Z_STREAM_END) {
                complete = true;
            } else if (status != Z_OK) {
                throw std::runtime_error("invalid gzip data or checksum: " + path.string());
            }
        }
        // Do not silently accept garbage or additional data after the member.
        if (stream.avail_in != 0 || source.peek() != std::char_traits<char>::eof()) {
            throw std::runtime_error("trailing data in gzip formatted checkpoint: " + path.string());
        }
    } catch (...) {
        inflateEnd(&stream);
        throw;
    }
    inflateEnd(&stream);
    return contents;
}

} // namespace

FormattedCheckpoint::FormattedCheckpoint(const std::filesystem::path& path) {
    std::error_code size_error;
    const auto file_size = std::filesystem::file_size(path, size_error);
    if (size_error) {
        throw std::runtime_error("cannot stat formatted checkpoint: " + path.string());
    }
    std::ifstream source(path, std::ios::binary);
    if (!source) {
        throw std::runtime_error("cannot open formatted checkpoint: " + path.string());
    }
    unsigned char signature[2]{};
    source.read(reinterpret_cast<char*>(signature), 2);
    const bool compressed = source.gcount() == 2 &&
                            signature[0] == 0x1f && signature[1] == 0x8b;
    source.clear();
    source.seekg(0);

    // Plain FCHK is streamed directly; gzip is decompressed into a buffer.
    // The array-size guard below must use the *uncompressed* size.
    std::string decoded;
    std::istringstream decoded_stream;
    std::istream* input_ptr = &source;
    std::uintmax_t logical_size = file_size;
    if (compressed) {
        source.close();
        decoded = decompress_gzip(path);
        logical_size = decoded.size();
        decoded_stream.str(std::move(decoded));
        input_ptr = &decoded_stream;
    }
    std::istream& input = *input_ptr;

    std::string line;
    if (!std::getline(input, line) || !std::getline(input, line)) {
        throw std::runtime_error("formatted checkpoint is missing its two-line header: " + path.string());
    }

    while (std::getline(input, line)) {
        if (line.size() < 44U) continue;
        const std::string label = trim_right(line.substr(0, 43));
        if (label.empty()) continue;

        std::istringstream meta(line.substr(43));
        char type = '\0';
        meta >> type;
        if (type != 'I' && type != 'R') continue;

        std::string marker;
        if (!(meta >> marker)) {
            throw std::runtime_error("malformed fchk record metadata: " + label);
        }

        Record record_value;
        if (marker == "N=") {
            long long count_raw = 0;
            if (!(meta >> count_raw) || count_raw < 0) {
                throw std::runtime_error("invalid array size in fchk: " + label);
            }
            const auto count_unsigned = static_cast<unsigned long long>(count_raw);
            if (count_unsigned > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max()) ||
                count_unsigned > static_cast<unsigned long long>(logical_size)) {
                throw std::runtime_error("unreasonable array size in fchk: " + label);
            }
            const auto count = static_cast<std::size_t>(count_unsigned);
            if (type == 'I') {
                record_value.kind = Kind::IntegerArray;
                record_value.integer_array.reserve(count);
                read_array_tokens(input, count, label, [&](const std::string& token) {
                    record_value.integer_array.push_back(parse_integer_token(token, label));
                });
            } else {
                record_value.kind = Kind::RealArray;
                record_value.real_array.reserve(count);
                read_array_tokens(input, count, label, [&](const std::string& token) {
                    record_value.real_array.push_back(parse_real_token(token, label));
                });
            }
        } else {
            if (type == 'I') {
                record_value.kind = Kind::IntegerScalar;
                record_value.integer_scalar = parse_integer_token(marker, label);
            } else {
                record_value.kind = Kind::RealScalar;
                record_value.real_scalar = parse_real_token(marker, label);
            }
        }

        const auto existing = records_.find(label);
        if (existing != records_.end()) {
            const auto& a = existing->second;
            const auto& b = record_value;
            bool identical = a.kind == b.kind;
            if (identical) {
                switch (a.kind) {
                case Kind::IntegerScalar:
                    identical = a.integer_scalar == b.integer_scalar;
                    break;
                case Kind::RealScalar:
                    identical = a.real_scalar == b.real_scalar;
                    break;
                case Kind::IntegerArray:
                    identical = a.integer_array == b.integer_array;
                    break;
                case Kind::RealArray:
                    identical = a.real_array == b.real_array;
                    break;
                }
            }
            if (!identical) {
                throw std::runtime_error("conflicting duplicate fchk record: " + label);
            }
            continue;
        }
        records_.emplace(label, std::move(record_value));
    }
}

const FormattedCheckpoint::Record& FormattedCheckpoint::record(const std::string& label) const {
    const auto it = records_.find(label);
    if (it == records_.end()) {
        throw std::runtime_error("missing fchk record: " + label);
    }
    return it->second;
}

bool FormattedCheckpoint::contains(const std::string& label) const {
    return records_.contains(label);
}

long long FormattedCheckpoint::integer(const std::string& label) const {
    const auto& r = record(label);
    if (r.kind != Kind::IntegerScalar) throw std::runtime_error("fchk record is not an integer scalar: " + label);
    return r.integer_scalar;
}

double FormattedCheckpoint::real(const std::string& label) const {
    const auto& r = record(label);
    if (r.kind != Kind::RealScalar) throw std::runtime_error("fchk record is not a real scalar: " + label);
    return r.real_scalar;
}

const std::vector<long long>& FormattedCheckpoint::integers(const std::string& label) const {
    const auto& r = record(label);
    if (r.kind != Kind::IntegerArray) throw std::runtime_error("fchk record is not an integer array: " + label);
    return r.integer_array;
}

const std::vector<double>& FormattedCheckpoint::reals(const std::string& label) const {
    const auto& r = record(label);
    if (r.kind != Kind::RealArray) throw std::runtime_error("fchk record is not a real array: " + label);
    return r.real_array;
}

} // namespace chem::kumo
