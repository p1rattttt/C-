#include<iostream>
#include<algorithm>
#include<cmath>
#include<vector>

template<typename T>
class Deque {
private:
    std::vector<T> container_;
    size_t l = 0;
    size_t r = 0;
    int shift = 0;
    void expand();
public:
    Deque() {
        container_.resize(4);
        l = 2;
        r = 2;
    }
    Deque(int len, T value = T()) {
        l = 0;
        r = len;
        container_.assign(len, value);
    }
    Deque(const Deque& x) {
        container_ = x.container_;
        l = x.l;
        r = x.r;
        shift = x.shift;
    }
    Deque& operator=(const Deque& x) {
        if (this == &x)
            return *this;
        container_ = x.container_;
        l = x.l;
        r = x.r;
        shift = x.shift;
        return *this;
    }

    size_t size() const;

    T& operator[](size_t ind);
    T operator[](size_t ind) const;

    T& at(size_t ind);
    T at(size_t ind) const;

    void push_back(T x);
    void push_front(T x);
    void pop_back();
    void pop_front();

    class iterator {
    protected:
        Deque<T>* nucl;
        int ind;
        int shift;
        size_t getIndex() const {
            return nucl->l + ind + nucl->shift - shift;
        }
    public:
        iterator() {
            nucl = nullptr;
            ind = 0;
            shift = 0;
        }

        iterator(Deque<T>* nucl, int ind, int shift) : nucl(nucl), ind(ind), shift(shift) {}

        iterator& operator++() {
            ++ind;
            return *this;
        }

        iterator operator++(int) {
            iterator ans = *this;
            ++ind;
            return ans;
        }

        iterator& operator--() {
            --ind;
            return *this;
        }

        iterator operator--(int) {
            iterator ans = *this;
            --ind;
            return ans;
        }

        iterator& operator+=(int x) {
            ind += x;
            return *this;
        }

        iterator& operator-=(int x) {
            ind -= x;
            return *this;
        }

        bool operator<(const iterator& x) const {
            return ind - shift < x.ind - x.shift;
        }

        bool operator==(const iterator& x) const {
            return !((*this) < x) && !(x < (*this));
        }

        bool operator!=(const iterator& x) const {
            return !((*this) == x);
        }

        bool operator>(const iterator& x) const {
            return x < (*this);
        }

        bool operator<=(const iterator& x) const {
            return !(x < (*this));
        }

        bool operator>=(const iterator& x) const {
            return !(x > (*this));
        }

        int operator-(const iterator& x) {
            return static_cast<int>(getIndex()) - x.getIndex();
        }

        iterator operator+(int x) {
            iterator ans = *this;
            ans += x;
            return ans;
        }

        iterator operator-(int x) {
            iterator ans = *this;
            ans -= x;
            return ans;
        }

        T& operator*() {
            return nucl->container_[getIndex()];
        }

        T* operator->() {
            return &nucl->container_[getIndex()];
        }
    };

    class const_iterator : public Deque<T>::iterator {
    public:
        const_iterator(Deque<T>* nucl, int ind, int shift) : Deque<T>::iterator(const_cast<Deque<T>*>(nucl), ind, shift) {}
        const T& operator*() const {
            return Deque<T>::iterator::nucl->container_[Deque<T>::iterator::getIndex()];
        }
        const T* operator->() const {
            return &Deque<T>::iterator::nucl->container_[Deque<T>::iterator::getIndex()];
        }
        
        const_iterator& operator++() {
            ++Deque<T>::iterator::ind;
            return *this;
        }

        const_iterator operator++(int) {
            iterator ans = *this;
            ++Deque<T>::iterator::ind;
            return ans;
        }

        const_iterator& operator--() {
            --Deque<T>::iterator::ind;
            return *this;
        }

        const_iterator operator--(int) {
            iterator ans = *this;
            --Deque<T>::iterator::ind;
            return ans;
        }

        const_iterator& operator+=(int x) {
            Deque<T>::iterator::ind += x;
            return *this;
        }

        const_iterator& operator-=(int x) {
            Deque<T>::iterator::ind -= x;
            return *this;
        }
		
        const_iterator operator+(int x) {
            const_iterator ans = *this;
            ans += x;
            return ans;
        }

        const_iterator operator-(int x) {
            const_iterator ans = *this;
            ans -= x;
            return ans;
        }
    };

    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;
    const_iterator cbegin() const;
    const_iterator cend() const;

    void insert(iterator it, T x);
    void erase(iterator it);
};

template<typename T>
void Deque<T>::expand() {
    size_t cursize = size();
    std::vector<T> copy = container_;
    container_.resize(3 * cursize);
    for (size_t i = l; i < r; ++i) {
        container_[i + static_cast<int>(cursize) - l] = copy[i];
    }
    r += static_cast<int>(cursize) - l;
    l += static_cast<int>(cursize) - l;
}

template<typename T>
size_t Deque<T>::size() const {
    return r - l;
}

template<typename T>
T& Deque<T>::operator[](size_t ind) {
    return container_[l + ind];
}

template<typename T>
T Deque<T>::operator[](size_t ind) const {
    return container_[l + ind];
}

template<typename T>
T& Deque<T>::at(size_t ind) {
    if (ind >= size())
        throw std::out_of_range("Segmention fault");
    return container_[l + ind];
}

template<typename T>
T Deque<T>::at(size_t ind) const {
    if (ind >= size())
        throw std::out_of_range("Segmention fault");
    return container_[l + ind];
}

template<typename T>
void Deque<T>::push_back(T x) {
    if (r + 1 >= container_.size())
        expand();
    container_[r++] = x;
}

template<typename T>
void Deque<T>::push_front(T x) {
    if (l == 0)
        expand();
    ++shift;
    container_[--l] = x;
}

template<typename T>
void Deque<T>::pop_back() {
    --r;
}

template<typename T>
void Deque<T>::pop_front() {
    --shift;
    ++l;
}

template<typename T>
typename Deque<T>::iterator Deque<T>::begin() {
    return Deque<T>::iterator(this, 0, shift);
}

template<typename T>
typename Deque<T>::iterator Deque<T>::end() {
    return Deque<T>::iterator(this, size(), shift);
}

template<typename T>
typename Deque<T>::const_iterator Deque<T>::cbegin() const {
    return Deque<T>::const_iterator(const_cast<Deque<T>*>(this), 0, shift);
}

template<typename T>
typename Deque<T>::const_iterator Deque<T>::cend() const {
    return Deque<T>::const_iterator(const_cast<Deque<T>*>(this), size(), shift);
}

template<typename T>
typename Deque<T>::const_iterator Deque<T>::begin() const {
    return cbegin();
}

template<typename T>
typename Deque<T>::const_iterator Deque<T>::end() const {
    return cend();
}

template<typename T>
void Deque<T>::insert(iterator it, T x) {
    push_back(x);
    for (auto i = --end(); i != it; --i) {
        std::swap(*i, *(i - 1));
    }
}

template<typename T>
void Deque<T>::erase(iterator it) {
    while (it != --end()) {
        std::swap(*it, *(it + 1));
        ++it;
    }
    pop_back();
}
