#ifndef FILTER_H_
#define FILTER_H_

#include <cstdint>
#include <string>
#include <vector>

namespace bench {

class Filter {
public:
    virtual std::string encode(const std::string& key, uint8_t* code_buf) = 0;
    virtual bool lookup(const std::string& key) = 0;
    virtual bool lookupRange(const std::string& left_key, const std::string& right_key) = 0;
    virtual uint64_t getMemoryUsage() = 0;
};

} // namespace bench

#endif // FILTER_H
