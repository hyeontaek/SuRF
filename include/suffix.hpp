#ifndef SUFFIX_H_
#define SUFFIX_H_

#include "bitvector.hpp"

#include <assert.h>

#include <algorithm>
#include <vector>
#include <unordered_set>

#include "config.hpp"
#include "hash.hpp"

namespace surf {

// Max suffix_len_ = 64 bits
// For kReal suffixes, if the stored key is not long enough to provide
// suffix_len_ suffix bits, its suffix field is cleared (i.e., all 0's)
// to indicate that there is no suffix info associated with the key.
class BitvectorSuffix : public Bitvector {
public:
    BitvectorSuffix() : type_(kNone), hash_suffix_len_(0), real_suffix_len_(0) {};

    BitvectorSuffix(const SuffixType type,
                    const level_t hash_suffix_len, const level_t real_suffix_len,
                    const std::vector<std::vector<word_t> >& bitvector_per_level,
                    const std::vector<position_t>& num_bits_per_level,
                    const level_t start_level = 0,
                    level_t end_level = 0/* non-inclusive */,
		    const std::vector<std::string>& intervals = std::vector<std::string>())
	: Bitvector(bitvector_per_level, num_bits_per_level, start_level, end_level) {
	assert((hash_suffix_len + real_suffix_len) <= kWordSize);
	type_ = type;
	hash_suffix_len_ = hash_suffix_len;
        real_suffix_len_ = real_suffix_len;
	intervals_ = intervals;
    }

    static word_t constructHashSuffix(const std::string& key, const level_t len) {
	word_t suffix = suffixHash(key);
	suffix <<= (kWordSize - len - kHashShift);
	suffix >>= (kWordSize - len);
	return suffix;
    }

    static word_t constructRealSuffix(const std::string& key,
				      const level_t level, const level_t len) {
	if (key.length() < level || ((key.length() - level) * 8) < len)
	    return 0;
	word_t suffix = 0;
	level_t num_complete_bytes = len / 8;
	if (num_complete_bytes > 0) {
	    suffix += (word_t)(label_t)key[level];
	    for (position_t i = 1; i < num_complete_bytes; i++) {
		suffix <<= 8;
		suffix += (word_t)(uint8_t)key[level + i];
	    }
	}
	level_t offset = len % 8;
	if (offset > 0) {
	    suffix <<= offset;
	    word_t remaining_bits = 0;
	    remaining_bits = (word_t)(uint8_t)key[level + num_complete_bytes];
	    remaining_bits >>= (8 - offset);
	    suffix += remaining_bits;
	}
	return suffix;
    }

    static word_t constructMixedSuffix(const std::string& key, const level_t hash_len,
				       const level_t real_level, const level_t real_len) {
        word_t hash_suffix = constructHashSuffix(key, hash_len);
        word_t real_suffix = constructRealSuffix(key, real_level, real_len);
        word_t suffix = hash_suffix;
        suffix <<= real_len;
        suffix |= real_suffix;
        return suffix;
    }

    static void findSuffixIntervals(std::vector<std::string>& intervals,
				    const std::vector<std::vector<const std::string*> >& keys_per_level,
				    const level_t real_len) {
	position_t count = 1 << static_cast<position_t>(real_len);

	// TODO: Sample suffixes to reduce sort time
	std::unordered_set<std::string> suffixes_set;
	for (level_t level = 0; level < keys_per_level.size(); level++) {
	    for (position_t i = 0; i < keys_per_level[level].size(); i++) {
		const std::string& key = *keys_per_level[level][i];
		if (key.size() < level)
		    suffixes_set.emplace();
		else
		    suffixes_set.emplace(key.c_str() + level, key.size() - level);
	    }
	}

	if (suffixes_set.size() == 0)
	    return;

	std::vector<std::string> suffixes(suffixes_set.begin(), suffixes_set.end());
	std::sort(suffixes.begin(), suffixes.end());

	intervals.clear();
	for (position_t i = 0; i < count - 1; i++) {
	    intervals.push_back(suffixes[suffixes.size() * i / count]);
	    //std::cout << intervals.back() << std::endl;
	}
	intervals.push_back(suffixes.back());
	//std::cout << intervals.back() << std::endl;
    }

    static word_t constructIntervalSuffix(const std::string& key,
					  const level_t real_level, const level_t real_len,
					  const std::vector<std::string>& intervals) {
	if (intervals.size() == 0)
	    return 0;

	// TODO: Avoid dynamic allocation
	std::string suffix;
	if (real_level < key.size())
	    suffix = key.substr(real_level);

	// XXX: Binary search can be replaced with another no-suffix SuRF traversal
	position_t left = 0;
	position_t right = intervals.size() - 1;
	position_t mid;
	while (left <= right) {
	    mid = (left + right) / 2;
	    int compare = intervals[mid].compare(suffix);
	    if (compare < 0)
		left = mid + 1;
	    else if (compare > 0)
		right = mid - 1;
	    else
		break;
	}

	// fixup -1 may not be necessary
	while (mid > 1 && intervals[mid - 1] > suffix) {
	    //std::cout << "fixup -1" << std::endl;
	    mid--;
	}

	while (mid + 1 < intervals.size() && suffix > intervals[mid]) {
	    //std::cout << "fixup +1" << std::endl;
	    mid++;
	}
	//std::cout << "fixup done" << std::endl;

	if (mid > 0)
	    assert(intervals[mid - 1] <= suffix);
	assert(mid < intervals.size());
	assert(suffix <= intervals[mid]);

	return static_cast<word_t>(mid);
    }

    static word_t constructSuffix(const SuffixType type, const std::string& key,
                                  const level_t hash_len,
                                  const level_t real_level, const level_t real_len) {
	switch (type) {
	case kHash:
	    return constructHashSuffix(key, hash_len);
	case kReal:
	    return constructRealSuffix(key, real_level, real_len);
        case kMixed:
            return constructMixedSuffix(key, hash_len, real_level, real_len);
        case kInterval:
	    // Must call constructIntervalSuffix() directly
	    assert(false);
            return 0;
	default:
	    return 0;
        }
    }

    static word_t extractHashSuffix(const word_t suffix, const level_t real_suffix_len) {
        return (suffix >> real_suffix_len);
    }

    static word_t extractRealSuffix(const word_t suffix, const level_t real_suffix_len) {
        word_t real_suffix_mask = 1;
        real_suffix_mask <<= real_suffix_len;
        real_suffix_mask--;
        return (suffix & real_suffix_mask);
    }

    SuffixType getType() const {
	return type_;
    }

    level_t getSuffixLen() const {
	return hash_suffix_len_ + real_suffix_len_;
    }

    level_t getHashSuffixLen() const {
	return hash_suffix_len_;
    }

    level_t getRealSuffixLen() const {
	return real_suffix_len_;
    }

    position_t serializedSize() const {
	position_t size = sizeof(num_bits_) + sizeof(type_)
            + sizeof(hash_suffix_len_) + sizeof(real_suffix_len_) + bitsSize();
	if (type_ == kInterval) {
	    for (size_t i = 0; i < intervals_.size(); i++) {
		size += sizeof(level_t);
		size += intervals_[i].size();
	    }
	}
	sizeAlign(size);
	return size;
    }

    position_t size() const {
	position_t size = (sizeof(BitvectorSuffix) + bitsSize());
	if (type_ == kInterval) {
	    for (size_t i = 0; i < intervals_.size(); i++) {
		// XXX: This may overestimate the memory use when short strings are inlined in std::string
		size += sizeof(intervals_[i]);
		size += intervals_[i].size();
	    }
	}
	return size;
    }

    word_t read(const position_t idx) const;
    word_t readReal(const position_t idx) const;
    bool checkEquality(const position_t idx, const std::string& key, const level_t level) const;

    // Compare stored suffix to querying suffix.
    // kReal suffix type only.
    int compare(const position_t idx, const std::string& key, const level_t level) const;

    int compareWithGreaterThanHint(const position_t idx, const std::string& key, const level_t level) const;
    int compareWithLessThanHint(const position_t idx, const std::string& key, const level_t level) const;

    void serialize(char*& dst) const {
	memcpy(dst, &num_bits_, sizeof(num_bits_));
	dst += sizeof(num_bits_);
	memcpy(dst, &type_, sizeof(type_));
	dst += sizeof(type_);
	memcpy(dst, &hash_suffix_len_, sizeof(hash_suffix_len_));
	dst += sizeof(hash_suffix_len_);
        memcpy(dst, &real_suffix_len_, sizeof(real_suffix_len_));
	dst += sizeof(real_suffix_len_);
	if (type_ != kNone) {
	    memcpy(dst, bits_, bitsSize());
	    dst += bitsSize();
	}
	if (type_ == kInterval) {
	    for (size_t i = 0; i < intervals_.size(); i++) {
		level_t len = static_cast<level_t>(intervals_[i].size());
		memcpy(dst, &len, sizeof(len));
		dst += sizeof(len);
		memcpy(dst, intervals_[i].c_str(), len);
		dst += len;
	    }
	}
	align(dst);
    }

    static BitvectorSuffix* deSerialize(char*& src) {
	BitvectorSuffix* sv = new BitvectorSuffix();
	memcpy(&(sv->num_bits_), src, sizeof(sv->num_bits_));
	src += sizeof(sv->num_bits_);
	memcpy(&(sv->type_), src, sizeof(sv->type_));
	src += sizeof(sv->type_);
	memcpy(&(sv->hash_suffix_len_), src, sizeof(sv->hash_suffix_len_));
	src += sizeof(sv->hash_suffix_len_);
        memcpy(&(sv->real_suffix_len_), src, sizeof(sv->real_suffix_len_));
	src += sizeof(sv->real_suffix_len_);
	if (sv->type_ != kNone) {
	    sv->bits_ = const_cast<word_t*>(reinterpret_cast<const word_t*>(src));
	    src += sv->bitsSize();
	}
	if (sv->type_ == kInterval) {
	    size_t count = 1 << static_cast<position_t>(sv->real_suffix_len_);
	    sv->intervals_.clear();
	    for (size_t i = 0; i < count; i++) {
		level_t len;
		memcpy(&len, src, sizeof(len));
		src += sizeof(len);
		sv->intervals_.emplace_back(src, len);
		src += len;
	    }
	}
	align(src);
	return sv;
    }

    void destroy() {
	if (type_ != kNone)
	    delete[] bits_;
    }

private:
    SuffixType type_;
    level_t hash_suffix_len_; // in bits
    level_t real_suffix_len_; // in bits
    std::vector<std::string> intervals_;
};

word_t BitvectorSuffix::read(const position_t idx) const {
    if (type_ == kNone) 
	return 0;

    level_t suffix_len = getSuffixLen();
    if (idx * suffix_len >= num_bits_) 
	return 0;

    position_t bit_pos = idx * suffix_len;
    position_t word_id = bit_pos / kWordSize;
    position_t offset = bit_pos & (kWordSize - 1);
    word_t ret_word = (bits_[word_id] << offset) >> (kWordSize - suffix_len);
    if (offset + suffix_len > kWordSize)
	ret_word += (bits_[word_id+1] >> (kWordSize - offset - suffix_len));
    return ret_word;
}

word_t BitvectorSuffix::readReal(const position_t idx) const {
    return extractRealSuffix(read(idx), real_suffix_len_);
}

bool BitvectorSuffix::checkEquality(const position_t idx, 
				    const std::string& key, const level_t level) const {
    if (type_ == kNone) 
	return true;
    if (idx * getSuffixLen() >= num_bits_) 
	return false;

    if (type_ == kInterval)
	return compareWithGreaterThanHint(idx, key, level) == kCouldBePositive;

    word_t stored_suffix = read(idx);
    if (type_ == kReal) {
	// if no suffix info for the stored key
	if (stored_suffix == 0) 
	    return true;
	// if the querying key is shorter than the stored key
	if (key.length() < level || ((key.length() - level) * 8) < real_suffix_len_) 
	    return false;
    }
    word_t querying_suffix 
	= constructSuffix(type_, key, hash_suffix_len_, level, real_suffix_len_);
    return (stored_suffix == querying_suffix);
}

// If no real suffix is stored for the key, compare returns 0.
// int BitvectorSuffix::compare(const position_t idx, 
// 			     const std::string& key, const level_t level) const {
//     if ((type_ == kNone) || (type_ == kHash) || (idx * getSuffixLen() >= num_bits_))
// 	return 0;
//     word_t stored_suffix = read(idx);
//     word_t querying_suffix = constructRealSuffix(key, level, real_suffix_len_);
//     if (type_ == kMixed)
//         stored_suffix = extractRealSuffix(stored_suffix, real_suffix_len_);

//     if (stored_suffix == 0) 
// 	return 0;
//     if (stored_suffix < querying_suffix) 
// 	return -1;
//     else if (stored_suffix == querying_suffix) 
// 	return 0;
//     else 
// 	return 1;
// }

int BitvectorSuffix::compare(const position_t idx, 
			     const std::string& key, const level_t level) const {
    if ((idx * getSuffixLen() >= num_bits_) || (type_ == kNone) || (type_ == kHash))
	return kCouldBePositive;

    if (type_ == kInterval)
	return compareWithGreaterThanHint(idx, key, level);

    word_t stored_suffix = read(idx);
    word_t querying_suffix = constructRealSuffix(key, level, real_suffix_len_);
    if (type_ == kMixed)
        stored_suffix = extractRealSuffix(stored_suffix, real_suffix_len_);

    if ((stored_suffix == 0) && (querying_suffix == 0))
	return kCouldBePositive;
    else if ((stored_suffix == 0) || (stored_suffix < querying_suffix))
	return -1;
    else if (stored_suffix == querying_suffix) 
	return kCouldBePositive;
    else 
	return 1;
}

int BitvectorSuffix::compareWithGreaterThanHint(const position_t idx, const std::string& key, const level_t level) const {
    if (type_ != kInterval)
	return compare(idx, key, level);

    word_t stored_suffix = read(idx);
    const std::string& lower_bound = intervals_[std::max(static_cast<int>(stored_suffix) - 1, 0)];

    std::string querying_suffix;
    if (level < key.size())
	querying_suffix = key.substr(level);

    int compare = lower_bound.compare(querying_suffix);
    if (compare <= 0) {
	const std::string& upper_bound = intervals_[stored_suffix];
	compare = upper_bound.compare(querying_suffix);
	if (compare >= 0)
	    compare = kCouldBePositive;
    }
    return compare;
}

int BitvectorSuffix::compareWithLessThanHint(const position_t idx, const std::string& key, const level_t level) const {
    if (type_ != kInterval)
	return compare(idx, key, level);

    word_t stored_suffix = read(idx);
    const std::string& upper_bound = intervals_[stored_suffix];

    std::string querying_suffix;
    if (level < key.size())
	querying_suffix = key.substr(level);

    int compare = upper_bound.compare(querying_suffix);
    if (compare >= 0) {
	const std::string& lower_bound = intervals_[std::max(static_cast<int>(stored_suffix) - 1, 0)];
	compare = lower_bound.compare(querying_suffix);
	if (compare <= 0)
	    compare = kCouldBePositive;
    }
    return compare;
}

} // namespace surf

#endif // SUFFIXVECTOR_H_
