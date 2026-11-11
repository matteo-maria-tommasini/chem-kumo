#pragma once

#include <array>
#include <filesystem>
#include <fstream>
#include <istream>
#include <streambuf>
#include <string>
#include <zlib.h>

namespace chem::kumo::detail {

/// True if the stream starts with the gzip signature 1f 8b. Compression is
/// detected from content, never from the file extension. The stream is
/// rewound to its beginning.
bool has_gzip_signature(std::istream& input);

/// Streaming decoder for a single-member gzip file, usable as the buffer of
/// a std::istream. Memory use is bounded (two 64 KiB buffers), and nothing is
/// written to a temporary file.
///
/// Integrity is checked strictly: the gzip trailer (CRC32 and uncompressed
/// size) must be reached and must match, and no data may follow the member.
/// Because iostream extraction turns streambuf failures into end-of-file,
/// errors are recorded and reported by throw_if_failed() and finish().
/// Error messages name the file and `description` (e.g. "Cube").
class GzipInputBuffer final : public std::streambuf {
public:
    GzipInputBuffer(const std::filesystem::path& path, std::string description);
    ~GzipInputBuffer() override;
    GzipInputBuffer(const GzipInputBuffer&) = delete;
    GzipInputBuffer& operator=(const GzipInputBuffer&) = delete;

    /// Throw std::runtime_error if a decoding error has occurred so far.
    void throw_if_failed() const;

    /// Consume the rest of the member and verify its trailer. Throws
    /// std::runtime_error for corrupt, truncated or trailing data.
    void finish();

protected:
    int_type underflow() override;

private:
    std::string message(const std::string& problem) const;

    std::filesystem::path path_;
    std::string description_;
    std::ifstream source_;
    z_stream stream_{};
    bool finished_ = false;
    std::string error_;
    std::array<unsigned char, 65536> input_{};
    std::array<char, 65536> output_{};
};

/// Decompress a whole gzip file into memory, with the checks of
/// GzipInputBuffer. Used where the uncompressed size must be known.
std::string read_gzip_file(const std::filesystem::path& path, const std::string& description);

} // namespace chem::kumo::detail
