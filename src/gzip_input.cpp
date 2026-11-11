#include "gzip_input.hpp"

#include <stdexcept>
#include <utility>

namespace chem::kumo::detail {

bool has_gzip_signature(std::istream& input) {
    unsigned char signature[2]{};
    input.read(reinterpret_cast<char*>(signature), 2);
    const bool gzip = input.gcount() == 2 && signature[0] == 0x1f && signature[1] == 0x8b;
    input.clear();
    input.seekg(0);
    return gzip;
}

GzipInputBuffer::GzipInputBuffer(const std::filesystem::path& path, std::string description)
    : path_(path), description_(std::move(description)), source_(path, std::ios::binary) {
    if (!source_) throw std::runtime_error(message("cannot open gzip"));
    // MAX_WBITS + 16: accept only the gzip wrapper, not raw zlib data.
    if (inflateInit2(&stream_, MAX_WBITS + 16) != Z_OK)
        throw std::runtime_error(message("cannot initialize gzip decoder for"));
    setg(output_.data(), output_.data(), output_.data());
}

GzipInputBuffer::~GzipInputBuffer() {
    inflateEnd(&stream_);
}

std::string GzipInputBuffer::message(const std::string& problem) const {
    return problem + " " + description_ + ": " + path_.string();
}

void GzipInputBuffer::throw_if_failed() const {
    if (!error_.empty()) throw std::runtime_error(error_);
}

void GzipInputBuffer::finish() {
    // Drain the remaining output: the trailer is verified only when inflate()
    // reaches the end of the member.
    while (underflow() != traits_type::eof())
        setg(eback(), egptr(), egptr());
    throw_if_failed();
    if (!finished_) throw std::runtime_error(message("truncated gzip"));
}

GzipInputBuffer::int_type GzipInputBuffer::underflow() {
    if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
    if (finished_ || !error_.empty()) return traits_type::eof();
    for (;;) {
        if (stream_.avail_in == 0) {
            source_.read(reinterpret_cast<char*>(input_.data()),
                         static_cast<std::streamsize>(input_.size()));
            const auto count = source_.gcount();
            if (count == 0) {
                error_ = message("truncated gzip");
                return traits_type::eof();
            }
            stream_.next_in = input_.data();
            stream_.avail_in = static_cast<uInt>(count);
        }
        stream_.next_out = reinterpret_cast<Bytef*>(output_.data());
        stream_.avail_out = static_cast<uInt>(output_.size());
        const int status = inflate(&stream_, Z_NO_FLUSH);
        const auto produced = output_.size() - stream_.avail_out;
        if (status == Z_STREAM_END) {
            finished_ = true;
            // Reject garbage or further members after the first one.
            if (stream_.avail_in != 0 || source_.peek() != std::char_traits<char>::eof())
                error_ = message("trailing data after gzip member in");
        } else if (status != Z_OK) {
            error_ = message("invalid gzip data or checksum in");
        }
        if (produced != 0) {
            setg(output_.data(), output_.data(), output_.data() + produced);
            return traits_type::to_int_type(*gptr());
        }
        if (finished_ || !error_.empty()) return traits_type::eof();
    }
}

std::string read_gzip_file(const std::filesystem::path& path, const std::string& description) {
    GzipInputBuffer buffer(path, description);
    std::string contents;
    std::array<char, 65536> chunk{};
    for (;;) {
        const auto count = buffer.sgetn(chunk.data(), static_cast<std::streamsize>(chunk.size()));
        if (count <= 0) break;
        contents.append(chunk.data(), static_cast<std::size_t>(count));
    }
    buffer.finish();
    return contents;
}

} // namespace chem::kumo::detail
