#include<cassert>
#include<iostream>
#include<cstring>
#include<cmath>

class String {
private:
    static const size_t MINIMAL_BUFFER_SIZE;
    static const double EXPANSION_COEFFICIENT;
    char* buffer_ = nullptr;
    size_t bufferSize_ = 0;
    size_t realSize_ = 0;
    void setNewSize_(size_t newSize);
public:
    const char& operator[](size_t index) const;
    char& operator[](size_t index);
    size_t length() const;
    void push_back(char c);
    void pop_back();
    char& front();
    const char& front() const;
    char& back();
    const char& back() const;
    size_t find(String sub) const;
    size_t rfind(String sub) const; 
    String substr(int start, int count) const;
    bool empty() const;
    void clear(); 
    void swap(String& a);
    String& operator=(const String& s);
    bool operator==(const String& s) const; 
    String& operator+=(const String &s);
    String& operator+=(char c);
    
    String() : buffer_(nullptr), bufferSize_(0), realSize_(0) {}

    String(size_t sz, char ch) : realSize_(sz) {
        setNewSize_(sz); 
        memset(buffer_, ch, realSize_);
    }

    String(const char* c) {
        setNewSize_(strlen(c));
        realSize_ = strlen(c);
        memcpy(buffer_, c, strlen(c));
    }
	
	String(char c) {
		setNewSize_(1);
		realSize_ = 1;
		buffer_[0] = c;
	}	

    ~String() {
        delete[] buffer_;
        bufferSize_ = 0;
        realSize_ = 0;
    }

    String(const String& s) {
        setNewSize_(s.length());
        realSize_ = s.length();
        memcpy(buffer_, s.buffer_, s.length());
    }
};

const size_t String::MINIMAL_BUFFER_SIZE = 4;
const double String::EXPANSION_COEFFICIENT = 2.0;

void String::setNewSize_(size_t newBufferSize) {
    newBufferSize = std::max(static_cast<size_t>(newBufferSize), MINIMAL_BUFFER_SIZE);
    if (newBufferSize == bufferSize_)
        return;
    char* newBuffer = new char[newBufferSize];
    if (bufferSize_) {
        std::memcpy(newBuffer, buffer_, sizeof(char) * realSize_);
        delete[] buffer_;
    }
    buffer_ = newBuffer;
    bufferSize_ = newBufferSize;
}

const char& String::operator[](size_t index) const {
    return buffer_[index];
}

char& String::operator[](size_t index) {
    return buffer_[index]; 
}

void String::push_back(char c) {
    if (bufferSize_ == realSize_)
        setNewSize_(std::ceil(bufferSize_ * EXPANSION_COEFFICIENT));
    buffer_[realSize_++] = c;
}

void String::pop_back() {
    --realSize_;
    if (realSize_ * EXPANSION_COEFFICIENT * EXPANSION_COEFFICIENT <= bufferSize_) 
        setNewSize_(std::floor(bufferSize_ / EXPANSION_COEFFICIENT));
}

size_t String::length() const {
    return realSize_;
}

void String::swap(String& a) {
    std::swap(buffer_, a.buffer_);
    std::swap(bufferSize_, a.bufferSize_);
    std::swap(realSize_, a.realSize_);
}

bool String::empty() const {
    return !length();
}

String& String::operator=(const String& s) {
    if (this == &s) 
        return *this;
    String copy = s;
    swap(copy);
    return *this; 
}

char& String::front() {
    return *buffer_;
}

const char& String::front() const {
    return *buffer_;
}

char& String::back() {
    return buffer_[realSize_ - 1];
}

const char& String::back() const {
    return buffer_[realSize_ - 1];
}

void String::clear() {
    setNewSize_(0);
    realSize_ = 0;
}

bool String::operator==(const String& s) const {
    if (length() != s.length()) return 0;
    for (size_t i = 0; i < s.length(); ++i) {
        if (buffer_[i] != s[i]) return 0;
    }
    return 1;
}

String& String::operator+=(const String& s) {
    if (length() + s.length() > bufferSize_)
        setNewSize_((length() + s.length()) * EXPANSION_COEFFICIENT);
    memcpy(buffer_ + length(), s.buffer_, s.length());
    realSize_ += s.length();
    return *this;
}

String& String::operator+=(char c) {
    push_back(c); 
    return *this;
}

String operator+(const String& a, const String& b) {
    String nw = a;
    nw += b;
    return nw;
}

String String::substr(int start, int count) const {
    String nw;
    for (size_t i = start; i < static_cast<size_t>(start + count); ++i)
        nw += buffer_[i];
    return nw;
}

size_t String::find(String s) const {
    size_t len = s.length();
    for (size_t i = 0; i + len <= length(); ++i) {
		if (!memcmp(buffer_ + i, s.buffer_, len)) {
			return i;
		}
    }
    return length();
}

size_t String::rfind(String s) const {
    size_t len = s.length();
    size_t mx = length();
    for (size_t i = 0; i + len <= length(); ++i) {
		if (!memcmp(buffer_ + i, s.buffer_, len)) {
			mx = i;
		}
    } 
    return mx;
}

std::ostream& operator<<(std::ostream& out, const String& s) {
    for (size_t i = 0; i < s.length(); ++i)
        out << s[i];
    return out;
}

std::istream& operator>>(std::istream& in, String& s) {
    s.clear();
    char c;
    in >> std::noskipws;
    while (in >> c && c != ' ' && c != '\n')
        s.push_back(c); 
    return in;
}

